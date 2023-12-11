#' Print the output of an \code{"otr_results"} object.
#'
#' @param x An \code{"otr_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print otr_results
#' @export
print.otr_results <- function(x, ...){

  res <- x$results

  threshold_names <- grep("^threshold", names(res), value=TRUE)

  for(t in 1:length(threshold_names)){
    threshold <- threshold_names[t]

    sub <- res[[threshold]]

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

#' Compare outcomes under different treatment rules
#'
#' @param res_rule1 An \code{"otr_results"} object for treatment rule 1
#' @param res_rule2 An \code{"otr_results"} object for treatment rule 2
#' @param threshold1 Threshold to use for treatment rule 1
#' @param threshold2 Threshold to use for treatment rule 2
#' @param rule1_comp Which effect for rule 1 to use in comparison ("treatment effect"/"te" or "subgroup effect"/"se")
#' @param rule2_comp Which effect for rule 2 to use in comparison ("treatment effect"/"te" or "subgroup effect"/"se")
#' @param ... Other arguments (not used)
#'
#' @import glue
#'
#' @returns dataframe containing expected value and variance for comparison
#' @export
compare.otr_results <- function(res_rule1, res_rule2, threshold1, threshold2, rule1_comp, rule2_comp, ...){
  # get name of threshold to pull results
  t_name1 <- glue("threshold =  ", threshold1)
  t_name2 <- glue("threshold =  ", threshold2)

  # get k fold assignments in each rule
  k_fold_assignments_rule1 <- res_rule1$results[[t_name1]]$decision_df
  k_fold_assignments_rule2 <- res_rule2$results[[t_name2]]$decision_df

  # get influence functions from each fold and bind into single matrix for each rule
  inf_fns_rule1 <- res_rule1$results[[t_name1]]$k_fold_results$influence_fns
  inf_fns_rule2 <- res_rule2$results[[t_name2]]$k_fold_results$influence_fns

  inf_fns_rule1 <- do.call(rbind, inf_fns_rule1)
  inf_fns_rule2 <- do.call(rbind, inf_fns_rule2)

  # reorder by participant ID number then drop id col
  inf_fns_rule1 <- cbind(id = k_fold_assignments_rule1$id, inf_fns_rule1)
  inf_fns_rule2 <- cbind(id = k_fold_assignments_rule2$id, inf_fns_rule2)

  id_order <- order(inf_fns_rule1[,"id"])

  inf_fns_rule1 <- inf_fns_rule1[id_order, ]
  inf_fns_rule2 <- inf_fns_rule2[id_order, ]

  inf_fns_rule1 <- inf_fns_rule1[,-c(1)]
  inf_fns_rule2 <- inf_fns_rule2[,-c(1)]

  # combine into single inf_fn_matrix
  inf_fn_matrix <- cbind(inf_fns_rule1, inf_fns_rule2)
  inf_fn_matrix <- apply(inf_fn_matrix, 2, as.numeric)

  # get covariance matrix
  cov_matrix <- stats::cov(inf_fn_matrix)

  # Pull AITPWs and make gradient for comparing rules
  if(rule1_comp == "treatment effect" | rule1_comp == "te"){
    aiptw_1 <- res_rule1$results[[t_name1]]$aggregated_results$treatment_effect
    gradient_1 <- c(0, 0, 0, 1, 0)
  } else if(rule1_comp == "subgroup effect" | rule1_comp == "se"){
    aiptw_1 <- res_rule1$results[[t_name1]]$aggregated_results$subgroup_effect
    gradient_1 <- c(0, 0, 0, 0, 1)
  } else{
    # error and exit
  }

  if(rule2_comp == "treatment effect" | rule2_comp == "te"){
    aiptw_2 <- res_rule2$results[[t_name2]]$aggregated_results$treatment_effect
    gradient_2 <- c(0, 0, 0, -1, 0)
  } else if(rule2_comp == "subgroup effect" | rule2_comp == "se"){
    aiptw_2 <- res_rule2$results[[t_name2]]$aggregated_results$subgroup_effect
    gradient_2 <- c(0, 0, 0, 0, -1)
  } else{
    # error and exit
  }

  gradient <- c(gradient_1, gradient_2)

  # calculate results
  exp_val_of_comparison <- aiptw_1 - aiptw_2
  var_of_comparison <- t(gradient) %*% cov_matrix %*% gradient

  return(data.frame(Z_list_1 = res_rule1$Z_list,
                    threshold_1 = threshold1,
                    Z_list_2 = res_rule2$Z_list,
                    threshold_2 = threshold2,
                    expected_val_of_comparison = exp_val_of_comparison,
                    var_of_comparison = var_of_comparison))
}
