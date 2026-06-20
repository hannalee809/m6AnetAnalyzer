#' Compute Wilcoxon Test (Mann-Whitney U test) for Differential Methylation
#'
#' This function computes a statistical test for differential methylation for
#' conditions with replicates. Specifically, it performs a Wilcoxon rank-sum
#' test per transcript to compare weighted modification ratios (WMR) between the
#' two conditions.
#'
#' @param condition1_list A list of data frames for condition 1. Each must contain
#'   `ensembl_transcript_id` and `weighted_mod_ratio`. This is the CONTROL group.
#' @param condition2_list A list of data frames for condition 2. Same required
#'   columns as `condition1_list`. This is the EXPERIMENTAL group.
#' @param condition1_name A string naming condition 1.
#' @param condition2_name A string naming condition 2.
#'
#' @return A tibble with one row per transcript containing:
#'   \itemize{
#'     \item `m6A transcript info`
#'     \item `median weighted mod ratio for condition 1`
#'     \item `median weighted mod ratio for condition 2`
#'     \item `median_diff` — median WMR( condition2 – condition1 )
#'     \item `p_value` — Wilcoxon rank-sum test p-value
#'   }
#'
#' @details Transcripts with all-zero WMR values across replicates are removed.
#'   Transcripts without data in both conditions return `NA` for all outputs.
#'
#' @importFrom dplyr bind_rows mutate group_by filter group_modify
#' @importFrom tibble tibble
#' @importFrom stats wilcox.test
#' @export

# Combine replicates and label them
run_wmr_differential_test <- function(condition1_list,
         condition2_list,
         condition1_name,
         condition2_name) {
  # Combine replicates and label them
  combine_reps <- function(df_list, condition_name) {
    bind_rows(lapply(seq_along(df_list), function(i) {
      df_list[[i]] %>% mutate(condition = condition_name, replicate = i)
    }))
  }

  # Combine all replicates
  all_wmr <- bind_rows(
    combine_reps(condition1_list, condition1_name),
    combine_reps(condition2_list, condition2_name)
  )

  # Filter out transcripts with all-zero weighted_mod_ratio
  filtered <- all_wmr %>%
    group_by(ensembl_transcript_id) %>%
    filter(any(weighted_mod_ratio > 0)) %>%
    dplyr::ungroup()

  # Run Wilcoxon per transcript and summarize
  results <- filtered %>%
    group_by(
      ensembl_transcript_name,
      ensembl_transcript_id,
      ensembl_gene_id,
      ensembl_gene_name,
      gene_biotype,
      transcript_length
    ) %>%
    group_modify( ~ {
      # Skip transcripts without both conditions
      if (length(unique(.x$condition)) < 2) {
        return(
          tibble(
            median_wmr_condition1 = NA,
            median_wmr_condition2 = NA,
            median_diff = NA,
            p_value = NA
          )
        )
      }

      median_c1 <-
        median(.x$weighted_mod_ratio[.x$condition == condition1_name])
      median_c2 <-
        median(.x$weighted_mod_ratio[.x$condition == condition2_name])
      wt <-
        wilcox.test(weighted_mod_ratio ~ condition,
                    data = .x,
                    exact = FALSE)

      tibble(
        median_wmr_condition1 = median_c1,
        median_wmr_condition2 = median_c2,
        median_diff = median_c2 - median_c1,
        p_value = wt$p.value
      )
    }) %>%
    dplyr::ungroup()

  return(results)
}
