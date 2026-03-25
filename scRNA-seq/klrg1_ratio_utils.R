# scRNA-seq/klrg1_ratio_utils.R
#
# Utilities for KLRG1-based cytotoxic/non-cytotoxic composition analysis.


compute_klrg1_cytotoxic_ratio <- function(
    metadata,
    patient_col,
    state_col,
    cytotoxic_label,
    klrg1_col = "KLRG1_expression",
    min_cells = 10,
    pseudocount = 1
) {
  required_cols <- c(patient_col, state_col, klrg1_col)
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  metadata %>%
    dplyr::mutate(
      ctl_status = dplyr::if_else(.data[[state_col]] == cytotoxic_label, "Cytotoxic", "Other")
    ) %>%
    dplyr::count(.data[[patient_col]], .data[[klrg1_col]], ctl_status, name = "ncells") %>%
    dplyr::rename(patient_id = .data[[patient_col]], KLRG1_expression = .data[[klrg1_col]]) %>%
    tidyr::complete(
      patient_id,
      KLRG1_expression = c("KLRG1+", "KLRG1-"),
      ctl_status = c("Cytotoxic", "Other"),
      fill = list(ncells = 0)
    ) %>%
    tidyr::pivot_wider(
      names_from = ctl_status,
      values_from = ncells,
      values_fill = 0
    ) %>%
    dplyr::filter((Cytotoxic + Other) >= min_cells) %>%
    dplyr::mutate(
      ratio = (Cytotoxic + pseudocount) / (Other + pseudocount),
      ratio_log2 = log2(ratio)
    )
}
