#' Print the output of an \code{"otr_results"} object.
#'
#' @param x An \code{"otr_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print otr_results
#' @export
print.otr_results <- function(x, ...){

  tmp <- data.frame(
    c(x$overall_results$aiptw_a_1, x$overall_results$aiptw_a_0, x$overall_results$treatment_effect),
    c(x$overall_results$se_aiptw_a_1, x$overall_results$se_aiptw_a_0, x$overall_results$se_treatment_effect),
    c(x$overall_results$aiptw_a_1 - 1.96*x$overall_results$se_aiptw_a_1,
      x$overall_results$aiptw_a_0 - 1.96*x$overall_results$se_aiptw_a_0,
      x$overall_results$treatment_effect - 1.96*x$overall_results$se_treatment_effect),
    c(x$overall_results$aiptw_a_1 + 1.96*x$overall_results$se_aiptw_a_1,
      x$overall_results$aiptw_a_0 + 1.96*x$overall_results$se_aiptw_a_0,
      x$overall_results$treatment_effect + 1.96*x$overall_results$se_treatment_effect)
  )

  row_names <- c("E[Y(d) | d(Z) = 1]", "E[Y(0) | d(Z) = 1]", "E[Y(0) - Y(1)]")
  col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")

  rownames(tmp) <- row_names
  colnames(tmp) <- col_names

  # Print header with dashed line
  cat(paste("                               Results Aggregated Across k = ", max(x$decision_df$k), " folds \n"))
  cat(paste(rep("-", 100), collapse = ""), "\n")
  cat(sprintf("%-25s%-20s%-20s%-20s%-20s\n", "", col_names[1], col_names[2], col_names[3], col_names[4]))
  cat(paste(rep("-", 100), collapse = ""), "\n")

  for(i in 1:nrow(tmp)){
    row_to_print <- tmp[i, ]

    # Adjust the widths as needed
    formatted_row <- sprintf("%-25s%-20s%-20s%-20s%-20s\n",
                             row.names(row_to_print),
                             round(row_to_print[1],4),
                             round(row_to_print[2],4),
                             round(row_to_print[3],4),
                             round(row_to_print[4],4))

    # Print the formatted row
    cat(paste(formatted_row))
  }

  decision_df <- x$decision_df
  prop_treated <- sum(decision_df$decision) / nrow(decision_df)

  cat(paste("\nProportion treated under OTR: ", round(prop_treated, 4)))

  Z_list <- x$CATE_models[[1]]$object$coefficients[-1]

  cat(paste("\nCovariates used in decision rule: ", paste(Z_list, collapse = ", ")))

  invisible(tmp)

}
