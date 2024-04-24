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

#' Helper function to remove unnecessary output from SuperLearner model to reduce output size
#'
#' Reference: https://www.r-bloggers.com/2019/12/need-to-save-rs-lm-or-glm-models-trim-the-fat/
#'
#' @param cm SuperLearner GLM or Earth model
#' @param earth flag for if model is SuperLearner Earth model, defaults to GLM
#'
#' @return model with unnecessary features removed
strip_glm <- function(cm, earth = FALSE) {
  cm$cvFitLibrary <- NULL
  cm$env <- NULL
  cm$y = c()
  cm$model = c()

  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
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

#' Helper function to apply strip_glm function to libraries in CATE model
#'
#' @param cate_model CATE model fit by SuperLearner
#'
#' @return cate_model of reduced output size
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

#' Helper function to apply strip_glm function to libraries in Nuisance model
#'
#' @param nuisance_model Nuisance model fit by SuperLearner
#'
#' @return Nuisance model of reduced output size
strip_nuisance <- function(nuisance_model){

  # eliminate cvmodels to save space in nuisance object
  nuisance_model$cvFitLibrary <- NULL

  # eliminate env to save space in nuisance object
  nuisance_model$env <- NULL

  for(i in 1:length(nuisance_model$fitLibrary)){

    if(class(nuisance_model$fitLibrary[[i]][[1]])[1] == "glm"){
      nuisance_model$fitLibrary[[i]][[1]] <- strip_glm(nuisance_model$fitLibrary[[i]][[1]])
    } else if(class(nuisance_model$fitLibrary[[i]])[1] == "SL.earth"){
      if(!is.null(nuisance_model$fitLibrary[[i]][[1]]$glm.list)){
        #check for presence of glm.list
        for(earth_glm in 1:length(nuisance_model$fitLibrary[[i]][[1]]$glm.list)){
          nuisance_model$fitLibrary[[i]][[1]]$glm.list[[earth_glm]] <- strip_glm(nuisance_model$fitLibrary[[i]][[1]]$glm.list[[earth_glm]], earth = TRUE)
        }

      }
    }
  }
  return(nuisance_model)
}

#' Print the output of a \code{"full_otr_results"} object.
#'
#' @param x A \code{"full_otr_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print full_otr_results
#' @export
print.full_otr_results <- function(x, ...){

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
        sub$aggregated_results$subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect),
      c(sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$se_treatment_effect),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 - 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect - 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 - 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect - 1.96*sub$aggregated_results$se_treatment_effect),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 + 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect + 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 + 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect + 1.96*sub$aggregated_results$se_treatment_effect)
    )

    row_names <- c("E[Y(d) | d(Z) = 1]",
                   "E[Y(0) | d(Z) = 1]",
                   "E[d(Z) = 1]",
                   "E[Y(d) - Y(0) | d(Z) = 1]",
                   "E[Y(1) - Y(d) | d(Z) = 0]",
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

#' Print the output of a \code{"full_otr_results"} object.
#'
#' @param x An \code{"otr_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print otr_results
#' @export
print.otr_results <- function(x, ...){

  res <- x

  threshold_names <- grep("^threshold", names(res), value=TRUE)

  for(t in 1:length(threshold_names)){
    threshold <- threshold_names[t]

    sub <- res[[threshold]]

    tmp <- data.frame(
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1,
        sub$aggregated_results$subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect),
      c(sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$se_treatment_effect),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 - 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect - 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 - 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect - 1.96*sub$aggregated_results$se_treatment_effect),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 + 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect + 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 + 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect + 1.96*sub$aggregated_results$se_treatment_effect)
    )

    row_names <- c("E[Y(d) | d(Z) = 1]",
                   "E[Y(0) | d(Z) = 1]",
                   "E[d(Z) = 1]",
                   "E[Y(d) - Y(0) | d(Z) = 1]",
                   "E[Y(1) - Y(d) | d(Z) = 0]",
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
#' @param res_rule1 An \code{"otr_results"} object or \code{"full_otr_results"} object for treatment rule 1
#' @param res_rule2 An \code{"otr_results"} object or \code{"full_otr_results"} object for treatment rule 2
#' @param threshold threshold to use for comparison of both rules (must appear in both otr_results objects)
#' @param rule1_comp Effect type for rule 1 to use in comparison ("treatment effect"/"te" or "subgroup effect"/"se")
#' @param rule2_comp Effect type for rule 2 to use in comparison ("treatment effect"/"te" or "subgroup effect"/"se")
#' @param ... Other arguments (not used)
#'
#' @returns dataframe containing expected value and variance for comparison
#' @export
compare.otr_results <- function(res_rule1, res_rule2, threshold, rule1_comp, rule2_comp, ...){

  # if just passed in otr_results object, put in list (so same format as full_otr_results)
  if(class(res_rule1) == "otr_results"){
    res_rule1 <- list(results = res_rule1)
  }

  if(class(res_rule2) == "otr_results"){
    res_rule2 <- list(results = res_rule2)
  }

  # get name of threshold to pull results
  t_name <- paste("threshold = ", threshold)

  # get k fold assignments in each rule
  k_fold_assignments_rule1 <- res_rule1$results[[t_name]]$decision_df
  k_fold_assignments_rule2 <- res_rule2$results[[t_name]]$decision_df

  # if either threshold wasn't found, error out
  if(is.null(k_fold_assignments_rule1) & is.null(k_fold_assignments_rule2)){
    stop(print("Thresholds for rule 1 and rule 2 not found"))
  } else if(is.null(k_fold_assignments_rule1)){
    stop(print("Threshold for rule 1 not found"))
  } else if(is.null(k_fold_assignments_rule2)){
    stop(print("Threshold for rule 2 not found"))
  }

  # get influence functions from each fold and bind into single matrix for each rule
  inf_fns_rule1 <- res_rule1$results[[t_name]]$k_fold_results$influence_fns
  inf_fns_rule2 <- res_rule2$results[[t_name]]$k_fold_results$influence_fns

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

  # Updated for subgroup effect in untreated
  if(rule1_comp == "treatment effect" | rule1_comp == "te"){
    aiptw_1 <- res_rule1$results[[t_name]]$aggregated_results$treatment_effect
    gradient_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)
  } else if(rule1_comp == "subgroup effect" | rule1_comp == "se"){
    aiptw_1 <- res_rule1$results[[t_name]]$aggregated_results$subgroup_effect
    gradient_1 <- c(0, 0, 0, 0, 0, 0, 1, 0, 0)
  } else if(rule1_comp == "subgroup effect untreated" | rule1_comp == "se_dZ0"){
    aiptw_1 <- res_rule1$results[[t_name]]$aggregated_results$subgroup_effect_dZ0
    gradient_1 <- c(0, 0, 0, 0, 0, 0, 0, 1, 0)
  } else{
    stop(print("Must enter treatment effect, te, subgroup effect, or se"))
  }

  if(rule2_comp == "treatment effect" | rule2_comp == "te"){
    aiptw_2 <- res_rule2$results[[t_name]]$aggregated_results$treatment_effect
    gradient_2 <- c(0, 0, 0, 0, 0, 0, 0, 0, -1)
  } else if(rule2_comp == "subgroup effect" | rule2_comp == "se"){
    aiptw_2 <- res_rule2$results[[t_name]]$aggregated_results$subgroup_effect
    gradient_2 <- c(0, 0, 0, 0, 0, 0, -1, 0, 0)
  } else if(rule1_comp == "subgroup effect untreated" | rule1_comp == "se_dZ0"){
    aiptw_2 <- res_rule2$results[[t_name]]$aggregated_results$subgroup_effect_dZ0
    gradient_2 <- c(0, 0, 0, 0, 0, 0, 0, -1, 0)
  } else {
    stop(print("Must enter treatment effect, te; subgroup effect, se; or subgroup effect untreated, se_dZ0"))
  }

  gradient <- c(gradient_1, gradient_2)

  # calculate results
  exp_val_of_comparison <- aiptw_1 - aiptw_2
  var_of_comparison <- t(gradient) %*% cov_matrix %*% gradient

  compare_rules <- data.frame(Z_list_1 = paste(res_rule1$results$Z_list, collapse=", "),
                              threshold = threshold,
                              rule1_comp = rule1_comp,
                              Z_list_2 = paste(res_rule2$results$Z_list, collapse=", "),
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
  } else if (x$rule1_comp == "se" | x$rule1_comp == "subgroup effect"){
    rule1_type <- "Subgroup Effect E[Y(d) - Y(0) | d(Z) = 1]"
  } else{
    rule1_type <- "Subgroup Effect in Untreated E[Y(1) - Y(d) | d(Z) = 0]"
  }

  threshold1 <- x$threshold

  rule1 <- paste(rule1_type, " for rule 1 at threshold = ", threshold1)

  # Make print statements for group 2
  if(x$rule2_comp == "te" | x$rule2_comp == "treatment effect"){
    rule2_type <- "Treatment Effect E[Y(d) - Y(0)]"
  } else if (x$rule2_comp == "se" | x$rule2_comp == "subgroup effect"){
    rule2_type <- "Subgroup Effect E[Y(d) - Y(0) | d(Z) = 1]"
  } else{
    rule2_type <- "Subgroup Effect in Untreated E[Y(1) - Y(d) | d(Z) = 0]"
  }

  threshold2 <- x$threshold

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

#' Helper function to average results across multiple seeds
#'
#' @param results_list list of results using same models & Z_list, different seeds
#' @param threshold threshold of results to use for comparison
#'
#' @return summary_df dataframe summarizing results across length(results_list) seeds
#' @export
#'
average_across_seeds <- function(results_list, threshold){

  # get name of threshold to use for comparison
  t_name <- paste("threshold = ", threshold)
  n <- length(results_list)
  combined_res <- data.frame()
  for(seed in 1:n){
    if(class(results_list[[seed]]) == "full_otr_results"){
      seed_res <- results_list[[seed]]$results[[t_name]]$aggregated_results
      if(is.null(seed_res)) stop(print(paste0("results_list item ", seed, " at threshold ", threshold, "not found")))
      Z_list <- results_list[[seed]]$results$Z_list

    } else if(class(results_list[[seed]]) == "otr_results"){
      seed_res <- results_list[[seed]][[t_name]]$aggregated_results
      if(is.null(seed_res)) stop(print(paste0("results_list item ", seed, " at threshold ", threshold, "not found")))
      Z_list <- results_list[[seed]]$Z_list

    } else {
      stop(print(paste0("results_list item ", seed, " not full_otr_results or otr_results object")))
    }

    # add results for given seed to dataframe
    combined_res <- rbind(combined_res, seed_res)

  }

  # average results across seeds
  average_results <- data.frame(
    EY_Ad_dZ1 = mean(combined_res$aiptw_EY_Ad_dZ1),
    se_EY_Ad_dZ1 = sqrt(mean((combined_res$se_aiptw_EY_Ad_dZ1^2))),
    EY_A0_dZ1 = mean(combined_res$aiptw_EY_A0_dZ1),
    se_EY_A0_dZ1 = sqrt(mean((combined_res$se_aiptw_EY_A0_dZ1^2))),
    E_dZ1 = mean(combined_res$E_dZ1),
    se_E_dZ1 = sqrt(mean((combined_res$se_E_dZ1^2))),
    subgroup_effect = mean(combined_res$subgroup_effect),
    se_subgroup_effect = sqrt(mean((combined_res$se_subgroup_effect^2))),
    subgroup_effect_dZ0 = mean(combined_res$subgroup_effect_dZ0),
    se_subgroup_effect_dZ0 = sqrt(mean((combined_res$se_subgroup_effect_dZ0^2))),
    treatment_effect = mean(combined_res$treatment_effect),
    se_treatment_effect = sqrt(mean((combined_res$se_treatment_effect^2)))
  )

  ci_df <- data.frame(
    lower_EY_Ad_dZ1 = average_results$EY_Ad_dZ1 - 1.96*average_results$se_EY_Ad_dZ1,
    upper_EY_Ad_dZ1 = average_results$EY_Ad_dZ1 + 1.96*average_results$se_EY_Ad_dZ1,
    lower_EY_A0_dZ1 = average_results$EY_A0_dZ1 - 1.96*average_results$se_EY_A0_dZ1,
    upper_EY_A0_dZ1 = average_results$EY_A0_dZ1 + 1.96*average_results$se_EY_A0_dZ1,
    lower_EY_E_dZ1 = average_results$E_dZ1 - 1.96*average_results$se_E_dZ1,
    upper_EY_E_dZ1 = average_results$E_dZ1 + 1.96*average_results$se_E_dZ1,
    lower_subgroup_effect = average_results$subgroup_effect - 1.96*average_results$se_subgroup_effect,
    upper_subgroup_effect = average_results$subgroup_effect + 1.96*average_results$se_subgroup_effect,
    lower_subgroup_effect_dZ0 = average_results$subgroup_effect_dZ0 - 1.96*average_results$se_subgroup_effect_dZ0,
    upper_subgroup_effect_dZ0 = average_results$subgroup_effect_dZ0 + 1.96*average_results$se_subgroup_effect_dZ0,
    lower_treatment_effect = average_results$treatment_effect - 1.96*average_results$se_treatment_effect,
    upper_treatment_effect = average_results$treatment_effect + 1.96*average_results$se_treatment_effect
  )

  results_table <- data.frame(
    value = c(average_results$EY_Ad_dZ1,
              average_results$EY_A0_dZ1,
              average_results$E_dZ1,
              average_results$subgroup_effect,
              average_results$subgroup_effect_dZ0,
              average_results$treatment_effect),
    se = c(average_results$se_EY_Ad_dZ1,
           average_results$se_EY_A0_dZ1,
           average_results$se_E_dZ1,
           average_results$se_subgroup_effect,
           average_results$se_subgroup_effect_dZ0,
           average_results$se_treatment_effect),
    lower_ci = c(ci_df$lower_EY_Ad_dZ1,
                 ci_df$lower_EY_A0_dZ1,
                 ci_df$lower_EY_E_dZ1,
                 ci_df$lower_subgroup_effect,
                 ci_df$lower_subgroup_effect_dZ0,
                 ci_df$lower_treatment_effect),
    upper_ci = c(ci_df$upper_EY_Ad_dZ1,
                 ci_df$upper_EY_A0_dZ1,
                 ci_df$upper_EY_E_dZ1,
                 ci_df$upper_subgroup_effect,
                 ci_df$upper_subgroup_effect_dZ0,
                 ci_df$upper_treatment_effect)
  )

  row_names <- c("E[Y(d) | d(Z) = 1]",
                 "E[Y(0) | d(Z) = 1]",
                 "E[d(Z) = 1]",
                 "E[Y(d) - Y(0) | d(Z) = 1]",
                 "E[Y(1) - Y(d) | d(Z) = 0]",
                 "E[Y(d) - Y(0)]")

  col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")

  rownames(results_table) <- row_names
  colnames(results_table) <- col_names

  average_results_list <- list(results_table = results_table,
                               Z_list = Z_list,
                               n_seeds = n,
                               threshold = threshold)

  class(average_results_list) <- "average_results"

  return(average_results_list)

}

#' Print the output of a \code{"average_results"} object.
#'
#' @param x An \code{"otr_results"} or \code{"full_otr_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print average_results
#' @export
print.average_results <- function(x, ...){
  Z_list <- x$Z_list
  n <- x$n_seeds
  results_table <- x$results_table

  row_names <- rownames(results_table)
  col_names <- colnames(results_table)

  # Print header with dashed line
  cat(paste("                                 Average results across n = ", n, " seeds \n"))
  cat(paste(rep("-", 105), collapse = ""), "\n")
  cat(sprintf("%-30s%-20s%-20s%-20s%-20s\n", "", col_names[1], col_names[2], col_names[3], col_names[4]))
  cat(paste(rep("-", 105), collapse = ""), "\n")

  for(i in 1:nrow(results_table)){
    row_to_print <- results_table[i, ]

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

  cat(paste("\nCovariates used in decision rule: ", paste(Z_list, collapse = ", "), "\n\n"))

}
