#' Method to average predictions over multiple SuperLearners
#' 
#' @param x Object of class \code{avgSuperLearner}
#' @param newdata Prediction \code{data.frame}
#' @param ... Other arguments (not used)
#' @return Vector of predictions on \code{newdata}
predict.avgSuperLearner <- function(x, newdata, ...){
  V <- length(x)
  pred_list <- lapply(x, predict, newdata = newdata)
  pred_list_sl <- lapply(pred_list, "[[", 1)
  avg_pred <- as.numeric(
    Reduce("+", pred_list_sl) / V
  )
  return(avg_pred)
}

# TO DO: Document
# helper function to eliminate unnecessary input in GLMs
# Reference: https://www.r-bloggers.com/2019/12/need-to-save-rs-lm-or-glm-models-trim-the-fat/
strip_glm <- function(cm, earth = FALSE) {
  cm$cvFitLibrary <- NULL
  cm$env <- NULL
  cm$y = c()
  cm$model = c()

  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()

  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()

  # handles earth package using offset from environment, environment dramatically increases size when saved to Rds
  if(earth == TRUE){
    attr(cm$terms,".Environment") = rlang::new_environment(data = list(offset = NULL), parent = baseenv())
  } else {
    attr(cm$terms,".Environment") = c()
  }

  attr(cm$formula,".Environment") = c()

  cm
}

# TO DO: Document
# helper function to apply strip_glm to any glm libraries in nuisance model
# else if to handle glm within earth
strip_cate <- function(cate_model){
  for(i in 1:length(cate_model$fitLibrary)){

    if(class(cate_model$fitLibrary[[i]][[1]])[1] == "glm"){
      cate_model$fitLibrary[[i]][[1]] <- strip_glm(cate_model$fitLibrary[[i]][[1]])
    } else if(class(cate_model$fitLibrary[[i]])[1] == "SL.earth"){
      if(!is.null(cate_model$fitLibrary[[i]][[1]]$glm.list)){
        #check for presence of glm.list
        for(earth_glm in 1:length(cate_model$fitLibrary[[i]][[1]]$glm.list)){
          cate_model$fitLibrary[[i]][[1]]$glm.list[[earth_glm]] <- strip_glm(cate_model$fitLibrary[[i]][[1]]$glm.list[[earth_glm]], earth = TRUE)
        }

      }
    }
  }
  return(cate_model)
}

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
#' @param threshold1 threshold to use for treatment rule 1
#' @param threshold2 threshold to use for treatment rule 2
#' @param rule1_comp Effect type for rule 1 to use in comparison ("treatment effect"/"te" or "subgroup effect"/"se")
#' @param rule2_comp Effect type for rule 2 to use in comparison ("treatment effect"/"te" or "subgroup effect"/"se")
#' @param ... Other arguments (not used)
#'
#' @returns dataframe containing expected value and variance for comparison
#' @export
compare.otr_results <- function(res_rule1, res_rule2, threshold1, threshold2, rule1_comp, rule2_comp, ...){
  # get name of threshold to pull results
  t_name1 <- paste("threshold = ", threshold1)
  t_name2 <- paste("threshold = ", threshold2)

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

  # If any of the influence functions are NA, drop those rows
  inf_fn_matrix <- inf_fn_matrix[!rowSums(is.na(inf_fn_matrix)),]

  # get covariance matrix
  cov_matrix <- stats::cov(inf_fn_matrix) / dim(inf_fn_matrix)[1]

  # Pull AIPTWs and make gradient for comparing rules
  if(rule1_comp == "treatment effect" | rule1_comp == "te"){
    aiptw_1 <- res_rule1$results[[t_name1]]$aggregated_results$treatment_effect
    gradient_1 <- c(0, 0, 0, 1, 0)
  } else if(rule1_comp == "subgroup effect" | rule1_comp == "se"){
    aiptw_1 <- res_rule1$results[[t_name1]]$aggregated_results$subgroup_effect
    gradient_1 <- c(0, 0, 0, 0, 1)
  } else{
    return(print("Must enter treatment effect, te, subgroup effect, or se"))
  }

  if(rule2_comp == "treatment effect" | rule2_comp == "te"){
    aiptw_2 <- res_rule2$results[[t_name2]]$aggregated_results$treatment_effect
    gradient_2 <- c(0, 0, 0, -1, 0)
  } else if(rule2_comp == "subgroup effect" | rule2_comp == "se"){
    aiptw_2 <- res_rule2$results[[t_name2]]$aggregated_results$subgroup_effect
    gradient_2 <- c(0, 0, 0, 0, -1)
  } else{
    return(print("Must enter treatment effect, te, subgroup effect, or se"))
  }

  gradient <- c(gradient_1, gradient_2)

  # calculate results
  exp_val_of_comparison <- aiptw_1 - aiptw_2
  var_of_comparison <- t(gradient) %*% cov_matrix %*% gradient

  compare_rules <- data.frame(Z_list_1 = paste(res_rule1$Z_list, collapse=", "),
                              threshold_1 = threshold1,
                              rule1_comp = rule1_comp,
                              Z_list_2 = paste(res_rule2$Z_list, collapse=", "),
                              threshold_2 = threshold2,
                              rule2_comp = rule2_comp,
                              expected_val_of_comparison = exp_val_of_comparison,
                              var_of_comparison = var_of_comparison)

  class(compare_rules) <- "otr_comparison"

  return(compare_rules)

}

#' Print the results of otr_comparison
#'
#' @param x An \code{"otr_comparison"} object
#' @param ... other arguments (not used)
#'
#' @method print otr_comparison
#' @export
print.otr_comparison <- function(x, ...){

  # Make print statements for group 1
  if(x$rule1_comp == "te" | x$rule1_comp == "treatment effect"){
    rule1_type <- "Treatment Effect E[Y(d) - Y(0)]"
  } else rule1_type <- "Subgroup Effect E[Y(d) - Y(0) | d(Z) = 1]"

  threshold1 <- x$threshold_1

  rule1 <- paste(rule1_type, " for rule 1 at threshold = ", threshold1)

  # Make print statements for group 2
  if(x$rule2_comp == "te" | x$rule2_comp == "treatment effect"){
    rule2_type <- "Treatment Effect E[Y(d) - Y(0)]"
  } else rule2_type <- "Subgroup Effect E[Y(d) - Y(0) | d(Z) = 1]"

  threshold2 <- x$threshold_2

  rule2 <- paste(rule2_type, " for rule 2 at threshold = ", threshold2)

  col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")

  cat(paste("\n", rule1, "\n", "vs \n", rule2, "\n"))
  cat(paste(rep("-", 105), collapse = ""), "\n")
  cat(sprintf("%-30s%-20s%-20s%-20s%-20s\n", "", col_names[1], col_names[2], col_names[3], col_names[4]))
  cat(paste(rep("-", 105), collapse = ""), "\n")
  cat(sprintf("%-30s%-20s%-20s%-20s%-20s\n",
              paste("Rule 1 - Rule 2"),
              round(x$expected_val_of_comparison, 4),
              round(sqrt(x$var_of_comparison), 4),
              round(x$expected_val_of_comparison - 1.96*sqrt(x$var_of_comparison), 4),
              round(x$expected_val_of_comparison + 1.96*sqrt(x$var_of_comparison), 4)))

  cat(paste("\n Rule 1: Z = ", x$Z_list_1))
  cat(paste("\n Rule 2: Z = ", x$Z_list_2))


}
