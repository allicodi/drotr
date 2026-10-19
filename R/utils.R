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
  
  threshold_names <- grep("^threshold", names(res), value = TRUE)
  
  for(t in 1:length(threshold_names)){
    threshold <- threshold_names[t]
    
    sub <- res[[threshold]]
    
    # Backwards compatibility for older objects
    new_fields <- c("E_Yd", "E_Y", "E_Yd_Y")
    all_new_fields <- c(new_fields, paste0("se_", new_fields))
    
    missing_fields <- setdiff(all_new_fields, names(sub$aggregated_results))
    
    # Fill absent fields with NA so data.frame construction works
    sub$aggregated_results[missing_fields] <- NA_real_
    
    tmp <- data.frame(
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1,
        sub$aggregated_results$subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect,
        sub$aggregated_results$compare_subgroup_effect,
        sub$aggregated_results$E_Yd,
        sub$aggregated_results$E_Y,
        sub$aggregated_results$E_Yd_Y),
      c(sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$se_treatment_effect,
        sub$aggregated_results$se_compare_subgroup_effect,
        sub$aggregated_results$se_E_Yd,
        sub$aggregated_results$se_E_Y,
        sub$aggregated_results$se_E_Yd_Y),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 - 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect - 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 - 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect - 1.96*sub$aggregated_results$se_treatment_effect,
        sub$aggregated_results$compare_subgroup_effect - 1.96*sub$aggregated_results$se_compare_subgroup_effect,
        sub$aggregated_results$E_Yd - 1.96*sub$aggregated_results$se_E_Yd,
        sub$aggregated_results$E_Y - 1.96*sub$aggregated_results$se_E_Y,
        sub$aggregated_results$E_Yd_Y - 1.96*sub$aggregated_results$se_E_Yd_Y),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 + 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect + 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 + 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect + 1.96*sub$aggregated_results$se_treatment_effect,
        sub$aggregated_results$compare_subgroup_effect + 1.96*sub$aggregated_results$se_compare_subgroup_effect,
        sub$aggregated_results$E_Yd + 1.96*sub$aggregated_results$se_E_Yd,
        sub$aggregated_results$E_Y + 1.96*sub$aggregated_results$se_E_Y,
        sub$aggregated_results$E_Yd_Y + 1.96*sub$aggregated_results$se_E_Yd_Y)
    )
    
    row_names <- c("E[Y(1) | d(Z) = 1]",
                   "E[Y(0) | d(Z) = 1]",
                   "E[d(Z) = 1]",
                   "E[Y(1) - Y(0) | d(Z) = 1]",
                   "E[Y(1) - Y(0) | d(Z) = 0]",
                   "E[Y(d) - Y(0)]",
                   "E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]",
                   "E[Y(d)]",
                   "E[Y]",
                   "E[Y(d) - Y]")
    
    col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")
    
    rownames(tmp) <- row_names
    colnames(tmp) <- col_names
    
    # Print header with dashed line
    cat(paste("                               Results for ", threshold,
              " Aggregated Across k = ", max(sub$k_non_na), " folds \n"))
    cat(paste(rep("-", 135), collapse = ""), "\n")
    cat(sprintf("%-60s%-20s%-20s%-20s%-20s\n",
                "", col_names[1], col_names[2], col_names[3], col_names[4]))
    cat(paste(rep("-", 135), collapse = ""), "\n")
    
    for(i in 1:nrow(tmp)){
      row_to_print <- tmp[i, ]
      
      vals <- as.character(round(unlist(row_to_print, use.names = FALSE), 4))
      
      # Print "--" only for fields absent from older objects
      if(i > 7){
        field <- new_fields[i - 7]
        
        if(field %in% missing_fields){
          vals[c(1, 3, 4)] <- "--"
        }
        
        if(paste0("se_", field) %in% missing_fields){
          vals[c(2, 3, 4)] <- "--"
        }
      }
      
      formatted_row <- sprintf("%-60s%-20s%-20s%-20s%-20s\n",
                               row.names(row_to_print),
                               vals[1],
                               vals[2],
                               vals[3],
                               vals[4])
      
      cat(formatted_row)
    }
    
    cat(paste("\nCovariates used in decision rule: ",
              paste(x$Z_list, collapse = ", "), "\n\n"))
  }
  
  invisible(tmp)
}

#' Print the output of an \code{"otr_results"} object.
#'
#' @param x An \code{"otr_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print otr_results
#' @export
print.otr_results <- function(x, ...){
  
  res <- x
  
  threshold_names <- grep("^threshold", names(res), value = TRUE)
  
  for(t in 1:length(threshold_names)){
    threshold <- threshold_names[t]
    
    sub <- res[[threshold]]
    
    # Backwards compatibility for older objects
    new_fields <- c("E_Yd", "E_Y", "E_Yd_Y")
    all_new_fields <- c(new_fields, paste0("se_", new_fields))
    
    missing_fields <- setdiff(all_new_fields, names(sub$aggregated_results))
    
    # Fill absent fields with NA so data.frame construction works
    sub$aggregated_results[missing_fields] <- NA_real_
    
    tmp <- data.frame(
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1,
        sub$aggregated_results$subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect,
        sub$aggregated_results$compare_subgroup_effect,
        sub$aggregated_results$E_Yd,
        sub$aggregated_results$E_Y,
        sub$aggregated_results$E_Yd_Y),
      c(sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$se_treatment_effect,
        sub$aggregated_results$se_compare_subgroup_effect,
        sub$aggregated_results$se_E_Yd,
        sub$aggregated_results$se_E_Y,
        sub$aggregated_results$se_E_Yd_Y),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 - 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 - 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect - 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 - 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect - 1.96*sub$aggregated_results$se_treatment_effect,
        sub$aggregated_results$compare_subgroup_effect - 1.96*sub$aggregated_results$se_compare_subgroup_effect,
        sub$aggregated_results$E_Yd - 1.96*sub$aggregated_results$se_E_Yd,
        sub$aggregated_results$E_Y - 1.96*sub$aggregated_results$se_E_Y,
        sub$aggregated_results$E_Yd_Y - 1.96*sub$aggregated_results$se_E_Yd_Y),
      c(sub$aggregated_results$aiptw_EY_Ad_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_Ad_dZ1,
        sub$aggregated_results$aiptw_EY_A0_dZ1 + 1.96*sub$aggregated_results$se_aiptw_EY_A0_dZ1,
        sub$aggregated_results$E_dZ1 + 1.96*sub$aggregated_results$se_E_dZ1,
        sub$aggregated_results$subgroup_effect + 1.96*sub$aggregated_results$se_subgroup_effect,
        sub$aggregated_results$subgroup_effect_dZ0 + 1.96*sub$aggregated_results$se_subgroup_effect_dZ0,
        sub$aggregated_results$treatment_effect + 1.96*sub$aggregated_results$se_treatment_effect,
        sub$aggregated_results$compare_subgroup_effect + 1.96*sub$aggregated_results$se_compare_subgroup_effect,
        sub$aggregated_results$E_Yd + 1.96*sub$aggregated_results$se_E_Yd,
        sub$aggregated_results$E_Y + 1.96*sub$aggregated_results$se_E_Y,
        sub$aggregated_results$E_Yd_Y + 1.96*sub$aggregated_results$se_E_Yd_Y)
    )
    
    row_names <- c("E[Y(1) | d(Z) = 1]",
                   "E[Y(0) | d(Z) = 1]",
                   "E[d(Z) = 1]",
                   "E[Y(1) - Y(0) | d(Z) = 1]",
                   "E[Y(1) - Y(0) | d(Z) = 0]",
                   "E[Y(d) - Y(0)]",
                   "E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]",
                   "E[Y(d)]",
                   "E[Y]",
                   "E[Y(d) - Y]")
    
    col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")
    
    rownames(tmp) <- row_names
    colnames(tmp) <- col_names
    
    # Print header with dashed line
    cat(paste("                               Results for ", threshold,
              " Aggregated Across k = ", max(sub$k_non_na), " folds \n"))
    cat(paste(rep("-", 135), collapse = ""), "\n")
    cat(sprintf("%-60s%-20s%-20s%-20s%-20s\n",
                "", col_names[1], col_names[2], col_names[3], col_names[4]))
    cat(paste(rep("-", 135), collapse = ""), "\n")
    
    for(i in 1:nrow(tmp)){
      row_to_print <- tmp[i, ]
      
      vals <- as.character(round(unlist(row_to_print, use.names = FALSE), 4))
      
      # Print "--" only for fields absent from older objects
      if(i > 7){
        field <- new_fields[i - 7]
        
        if(field %in% missing_fields){
          vals[c(1, 3, 4)] <- "--"
        }
        
        if(paste0("se_", field) %in% missing_fields){
          vals[c(2, 3, 4)] <- "--"
        }
      }
      
      formatted_row <- sprintf("%-60s%-20s%-20s%-20s%-20s\n",
                               row.names(row_to_print),
                               vals[1],
                               vals[2],
                               vals[3],
                               vals[4])
      
      cat(formatted_row)
    }
    
    cat(paste("\nCovariates used in decision rule: ",
              paste(x$Z_list, collapse = ", "), "\n\n"))
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
    gradient_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
  } else if(rule1_comp == "subgroup effect" | rule1_comp == "se"){
    aiptw_1 <- res_rule1$results[[t_name]]$aggregated_results$subgroup_effect
    gradient_1 <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
  } else if(rule1_comp == "subgroup effect untreated" | rule1_comp == "se_dZ0"){
    aiptw_1 <- res_rule1$results[[t_name]]$aggregated_results$subgroup_effect_dZ0
    gradient_1 <- c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
  } else{
    stop(print("Must enter treatment effect, te, subgroup effect, or se"))
  }

  if(rule2_comp == "treatment effect" | rule2_comp == "te"){
    aiptw_2 <- res_rule2$results[[t_name]]$aggregated_results$treatment_effect
    gradient_2 <- c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0)
  } else if(rule2_comp == "subgroup effect" | rule2_comp == "se"){
    aiptw_2 <- res_rule2$results[[t_name]]$aggregated_results$subgroup_effect
    gradient_2 <- c(0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0)
  } else if(rule2_comp == "subgroup effect untreated" | rule2_comp == "se_dZ0"){
    aiptw_2 <- res_rule2$results[[t_name]]$aggregated_results$subgroup_effect_dZ0
    gradient_2 <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0)
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
  
  # Backwards compatibility for older objects
  new_fields <- c("E_Yd", "E_Y", "E_Yd_Y")
  all_new_fields <- c(new_fields, paste0("se_", new_fields))
  
  missing_fields <- setdiff(all_new_fields, names(combined_res))
  
  # Fill absent fields with NA for calculations
  combined_res[missing_fields] <- NA_real_
  
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
    se_treatment_effect = sqrt(mean((combined_res$se_treatment_effect^2))),
    compare_subgroup_effect = mean(combined_res$compare_subgroup_effect),
    se_compare_subgroup_effect = sqrt(mean((combined_res$se_compare_subgroup_effect^2))),
    E_Yd = mean(combined_res$E_Yd),
    se_E_Yd = sqrt(mean((combined_res$se_E_Yd^2))),
    E_Y = mean(combined_res$E_Y),
    se_E_Y = sqrt(mean((combined_res$se_E_Y^2))),
    E_Yd_Y = mean(combined_res$E_Yd_Y),
    se_E_Yd_Y = sqrt(mean((combined_res$se_E_Yd_Y^2)))
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
    upper_treatment_effect = average_results$treatment_effect + 1.96*average_results$se_treatment_effect,
    lower_compare_subgroup_effect = average_results$compare_subgroup_effect - 1.96*average_results$se_compare_subgroup_effect,
    upper_compare_subgroup_effect = average_results$compare_subgroup_effect + 1.96*average_results$se_compare_subgroup_effect,
    lower_E_Yd = average_results$E_Yd - 1.96*average_results$se_E_Yd,
    upper_E_Yd = average_results$E_Yd + 1.96*average_results$se_E_Yd,
    lower_E_Y = average_results$E_Y - 1.96*average_results$se_E_Y,
    upper_E_Y = average_results$E_Y + 1.96*average_results$se_E_Y,
    lower_E_Yd_Y = average_results$E_Yd_Y - 1.96*average_results$se_E_Yd_Y,
    upper_E_Yd_Y = average_results$E_Yd_Y + 1.96*average_results$se_E_Yd_Y
  )
  
  results_table <- data.frame(
    value = c(average_results$EY_Ad_dZ1,
              average_results$EY_A0_dZ1,
              average_results$E_dZ1,
              average_results$subgroup_effect,
              average_results$subgroup_effect_dZ0,
              average_results$treatment_effect,
              average_results$compare_subgroup_effect,
              average_results$E_Yd,
              average_results$E_Y,
              average_results$E_Yd_Y),
    se = c(average_results$se_EY_Ad_dZ1,
           average_results$se_EY_A0_dZ1,
           average_results$se_E_dZ1,
           average_results$se_subgroup_effect,
           average_results$se_subgroup_effect_dZ0,
           average_results$se_treatment_effect,
           average_results$se_compare_subgroup_effect,
           average_results$se_E_Yd,
           average_results$se_E_Y,
           average_results$se_E_Yd_Y),
    lower_ci = c(ci_df$lower_EY_Ad_dZ1,
                 ci_df$lower_EY_A0_dZ1,
                 ci_df$lower_EY_E_dZ1,
                 ci_df$lower_subgroup_effect,
                 ci_df$lower_subgroup_effect_dZ0,
                 ci_df$lower_treatment_effect,
                 ci_df$lower_compare_subgroup_effect,
                 ci_df$lower_E_Yd,
                 ci_df$lower_E_Y,
                 ci_df$lower_E_Yd_Y),
    upper_ci = c(ci_df$upper_EY_Ad_dZ1,
                 ci_df$upper_EY_A0_dZ1,
                 ci_df$upper_EY_E_dZ1,
                 ci_df$upper_subgroup_effect,
                 ci_df$upper_subgroup_effect_dZ0,
                 ci_df$upper_treatment_effect,
                 ci_df$upper_compare_subgroup_effect,
                 ci_df$upper_E_Yd,
                 ci_df$upper_E_Y,
                 ci_df$upper_E_Yd_Y)
  )
  
  row_names <- c("E[Y(1) | d(Z) = 1]",
                 "E[Y(0) | d(Z) = 1]",
                 "E[d(Z) = 1]",
                 "E[Y(1) - Y(0) | d(Z) = 1]",
                 "E[Y(1) - Y(0) | d(Z) = 0]",
                 "E[Y(d) - Y(0)]",
                 "E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]",
                 "E[Y(d)]",
                 "E[Y]",
                 "E[Y(d) - Y]")
  
  col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")
  
  rownames(results_table) <- row_names
  colnames(results_table) <- col_names
  
  average_results_list <- list(
    results_table = results_table,
    Z_list = Z_list,
    n_seeds = n,
    threshold = threshold,
    missing_fields = missing_fields
  )
  
  class(average_results_list) <- "average_results"
  
  return(average_results_list)
}

#' Print the output of a \code{"average_results"} object.
#'
#' @param x An \code{"average_results"} object.
#' @param ... Other arguments (not used)
#'
#' @method print average_results
#' @export
print.average_results <- function(x, ...){
  
  Z_list <- x$Z_list
  n <- x$n_seeds
  results_table <- x$results_table
  threshold <- x$threshold
  
  # Backwards compatibility for older objects
  new_rows <- c("E[Y(d)]", "E[Y]", "E[Y(d) - Y]")
  new_fields <- c("E_Yd", "E_Y", "E_Yd_Y")
  
  # Rows absent from previously saved objects
  missing_rows <- setdiff(new_rows, rownames(results_table))
  
  if(length(missing_rows) > 0){
    missing_df <- as.data.frame(
      matrix(NA_real_,
             nrow = length(missing_rows),
             ncol = ncol(results_table),
             dimnames = list(missing_rows, colnames(results_table)))
    )
    
    results_table <- rbind(results_table, missing_df)
  }
  
  # Identify fields missing when averaging old seed objects
  missing_fields <- x$missing_fields
  
  row_names <- rownames(results_table)
  col_names <- colnames(results_table)
  
  # Print header with dashed line
  cat(paste("                                        Average results across n = ",
            n, " seeds for threshold ", threshold, "\n"))
  cat(paste(rep("-", 135), collapse = ""), "\n")
  cat(sprintf("%-60s%-20s%-20s%-20s%-20s\n",
              "", col_names[1], col_names[2], col_names[3], col_names[4]))
  cat(paste(rep("-", 135), collapse = ""), "\n")
  
  for(i in 1:nrow(results_table)){
    row_to_print <- results_table[i, ]
    
    vals <- as.character(round(unlist(row_to_print, use.names = FALSE), 4))
    
    # Print "--" for absent fields, preserving genuine NA values
    if(row_names[i] %in% new_rows){
      
      j <- match(row_names[i], new_rows)
      field <- new_fields[j]
      
      # Entire row absent from previously saved object
      if(row_names[i] %in% missing_rows){
        vals[] <- "--"
      } else {
        
        # Estimate absent from original seed objects
        if(field %in% missing_fields){
          vals[c(1, 3, 4)] <- "--"
        }
        
        # Standard error absent from original seed objects
        if(paste0("se_", field) %in% missing_fields){
          vals[c(2, 3, 4)] <- "--"
        }
      }
    }
    
    formatted_row <- sprintf("%-60s%-20s%-20s%-20s%-20s\n",
                             row_names[i],
                             vals[1],
                             vals[2],
                             vals[3],
                             vals[4])
    
    cat(formatted_row)
  }
  
  cat(paste("\nCovariates used in decision rule: ",
            paste(Z_list, collapse = ", "), "\n\n"))
}
