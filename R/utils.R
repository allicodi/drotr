#' Print the output of an \code{"otr_results"} object.
#'
#' @param x An \code{"otr_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print otr_results
#' @export
print.otr_results <- function(x, ...){

  threshold_names <- grep("^threshold", names(x), value=TRUE)

  for(t in 1:length(threshold_names)){
    threshold <- threshold_names[t]

    sub <- x[[threshold]]

    tmp <- data.frame(
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1,
        sub$aggregated_results$subgroup_effect,
        sub$aggregated_results$treatment_effect),
      c(sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$se_treatment_effect),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 - 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect - 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$treatment_effect - 1.96*sub$aggregated_results$se_treatment_effect),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 + 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect + 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$treatment_effect + 1.96*sub$aggregated_results$se_treatment_effect)
    )

    row_names <- c("E[Y(d) | d(Z) = 1]",
                   "E[Y(0) | d(Z) = 1]",
                   "E[d(Z) = 1]",
                   "E[Y(d) - Y(0) | d(Z) = 1]",
                   "E[Y(d) - Y(0)]")

    col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")

    rownames(tmp) <- row_names
    colnames(tmp) <- col_names

    # Print header with dashed line
    cat(paste("                      Results for ", threshold, " Aggregated Across k = ", max(sub$decision_df$k), " folds \n"))
    cat(paste(rep("-", 105), collapse = ""), "\n")
    cat(sprintf("%-30s%-20s%-20s%-20s%-20s\n", "", col_names[1], col_names[2], col_names[3], col_names[4]))
    cat(paste(rep("-", 105), collapse = ""), "\n")

    for(i in 1:nrow(tmp)){
      row_to_print <- tmp[i, ]

      # Adjust the widths as needed
      formatted_row <- sprintf("%-30s%-20s%-20s%-20s%-20s\n",
                               row.names(row_to_print),
                               round(row_to_print[1],4),
                               round(row_to_print[2],4),
                               round(row_to_print[3],4),
                               round(row_to_print[4],4))

      # Print the formatted row
      cat(paste(formatted_row))
    }

    cat(paste("\nCovariates used in decision rule: ", paste(x$Z_list, collapse = ", "), "\n\n"))

  }

  invisible(tmp)

}
