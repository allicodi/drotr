#' Function to apply existing OTR(s) to external data
#' 
#' @param df dataframe containing external dataset to apply rule(s) to; must contain Z_list variables that are same as pre-trained rule(s)
#' @param CATE_models list of CATE model(s) to apply to external dataset
#' @param Y_name name of outcome variable. Outcome variable in rule being applied should be the same as outcome variable in new data. 
#' @param A_name name of treatment variable. Treatment variable in rule being applied should be the same as the treatment variable in the new data. 
#' @param W_list character vector containing names of covariates in the dataframe used for nuisance models. 
#' @param Z_list character vector containing names of variables in df used to fit CATE model (variables used in treatment rule; must be same names as used in pre-fit CATE model(s))
#' @param id_name name of participant ID variable if present in data
#' @param nuisance_models list of objects of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models (only include if using pre-fit nuisance models)
#' @param sl.library.outcome character vector of SuperLearner libraries to use to fit the outcome models (if fitting nuisance internally)
#' @param sl.library.treatment character vector of SuperLearner libraries to use to fit the treatment models  (if fitting nuisance internally)
#' @param sl.library.missingness character vector of SuperLearner libraries to use to fit the missingness models  (if fitting nuisance internally)
#' @param threshold character vector of decision thresholds for CATE to determine OTR. Values should be positive if `Y_name` is desirable outcome, negative if `Y_name` is undesirable outcome. If threshold is 0, use +0 for desirable, -0 for undesirable.
#' @param k_folds integer number of folds to use for nuisance model cross-validation (must specify if fitting nuisance models here)
#' @param ps_trunc_level numeric level below which propensity scores will be truncated (to avoid errors in computing AIPTW)
#' @param outcome_type outcome_type specifying continuous (outcome_type = "gaussian") or binary (outcome_type = "binomial") outcome Y (if not providing pre-fit nuisance models)
#' @param truncate_CATE logical to indicate if large CATE predictions should be truncated at -1 and 1 (default = TRUE)
#'
#' @export
#'
#' @returns
#' \describe{
#'  \item{\code{results}}{list of `Results` objects for each threshold. See description of `Results` object in `compute_estimates`}
#'  \item{\code{nuisance_models}}{list of `Nuisance` objects containing outcome, treatment, and missingness models used in each fold}
#'  \item{\code{CATE_models}}{CATE model used in each fold}
#'  \item{\code{Z_list}}{character vector containing names of variables in df used to fit CATE model (variables used in treatment rule)}
#'  }
apply_OTR <- function(df,
                       CATE_models,
                       Y_name,
                       A_name,
                       W_list,
                       Z_list,
                       id_name = NULL,
                       threshold = c("0.05"),
                       nuisance_models = NULL,
                       sl.library.outcome = NULL,
                       sl.library.treatment = NULL,
                       sl.library.missingness = NULL,
                       k_folds = 5,
                       ps_trunc_level = 0.01,
                       outcome_type = "gaussian",
                       truncate_CATE = "TRUE"){
  
  # --------------------------------------------------------------------------
  # 0 - Prep dataset
  # --------------------------------------------------------------------------
  
  # if id_name not specified, make index id
  if(is.null(id_name)){
    df$id <- 1:nrow(df)
  } else {
    df$id <- df[[id_name]]
  }
  
  # --------------------------------------------------------------------------
  # 1 - Predict from each CATE model on full dataset
  # --------------------------------------------------------------------------
  CATE_preds <- drotr::predict_CATE_external(df, CATE_models, Z_list, truncate_CATE)
  
  # --------------------------------------------------------------------------
  # 2 - Make treatment decisions and compute treatment effects
  # --------------------------------------------------------------------------
  results <- drotr::compute_estimates_external(df, 
                                               Y_name, 
                                               A_name, 
                                               W_list, 
                                               Z_list, 
                                               CATE_preds, 
                                               threshold,
                                               nuisance_models,
                                               sl.library.outcome,
                                               sl.library.treatment,
                                               sl.library.missingness,
                                               k_folds,
                                               ps_trunc_level,
                                               outcome_type, 
                                               id_name)
  
  results <- list(results)
  names(results) <- "results"
  
  results$nuisance_models <- nuisance_models
  results$CATE_models <- CATE_models
  results$Z_list <- Z_list
  results$results$Z_list <- Z_list # adding in two places so full object print or just results subsection prints (so could opt to save partial results)
  
  class(results$results) <- "otr_results"
  class(results) <- "full_otr_results"
  
  return(results)
  
}

#' Function to predict CATEs & get treatment decisions on external dataset
#' 
#' @param df dataframe containing external dataset to apply rule(s) to; must contain Z_list variables that are same as pre-trained rule(s)
#' @param CATE_models list of CATE model(s) to apply to external dataset
#' @param Z_list character vector containing names of variables in df used to fit CATE model (variables used in treatment rule; must be same names as used in pre-fit CATE model(s))
#' @param truncate_CATE logical to indicate if large CATE predictions should be truncated at -1 and 1 (default = TRUE)
#' 
#' @export
#' 
#' @returns dataframe with columns for id, CATEs for each model, and average CATE for rule
predict_CATE_external <- function(df, CATE_models, Z_list, truncate_CATE = TRUE){
  
  Z <- df[, Z_list, drop = FALSE]
  
  # Get predicted CATEs for each model 
  CATE_preds <- lapply(CATE_models, function(CATE_model, truncate_CATE){
    # get prediction from CATE model 
    CATE_pred <- stats::predict(CATE_model, Z, type = 'response')
    
    # if truncate CATE is true, truncate at -1 and 1
    if(truncate_CATE == TRUE){
      trunc_flag <- ifelse(CATE_pred < -1, 1, 0)
      CATE_pred <- ifelse(CATE_pred < -1, -1, CATE_pred)
      
      trunc_flag <- ifelse(CATE_pred > 1, 1, trunc_flag)
      CATE_pred <- ifelse(CATE_pred > 1, 1, CATE_pred)
    } else {
      trunc_flag <- rep(0, length(CATE_pred))
    }
    
    return(data.frame(CATE_pred = CATE_pred, trunc_flag = trunc_flag))
  }, truncate_CATE = truncate_CATE)
  
  #avg_CATE <- rowMeans(data.frame(CATE_preds), na.rm = TRUE)
  CATE_preds <- data.frame(id = df$id, data.frame(CATE_preds))
  
  # colnames id, model_x, trunc_flag_x
  colnames(CATE_preds) <- c("id", 
                            unlist(lapply(1:length(CATE_models), function(i) {
                              c(paste0("model_", i), paste0("trunc_flag_", i))
                              })))
  
  model_cols <- grep("^model_", colnames(CATE_preds), value = TRUE)
  CATE_preds$avg_CATE <- rowMeans(CATE_preds[, model_cols], na.rm = TRUE)
  
  return(CATE_preds)

}

#' Function to compute estimates on external dataset
#' 
#' @param df dataframe containing external dataset to apply rule(s) to; must contain Z_list variables that are same as pre-trained rule(s)
#' @param CATE_preds dataframe with CATE predictions for each CATE model & truncation flags for extreme predictions
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param W_list character vector containing names of covariates in the dataframe to be used for fitting nuisance models
#' @param Z_list character vector containing names of variables in df used to fit CATE model (variables used in treatment rule; must be same names as used in pre-fit CATE model(s))
#' @param threshold character vector of decision thresholds for CATE to determine OTR. Values should be positive if `Y_name` is desirable outcome, negative if `Y_name` is undesirable outcome. If threshold is 0, use +0 for desirable, -0 for undesirable.
#' @param te_method method for estimating treatment effects, default aiptw
#' @param nuisance_models list of objects of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models 
#' 
#' @export
#' 
#' @returns \describe{
#' List of "Results" objects for each threshold. Each results object contains:
#' \item{\code{aggregated_results}}{Results averaged over each nuisance model}
#' \item{\code{k_fold_results}}{Results for each nuisance model}
#' \item{\code{decision_df}}{Dataframe with CATE prediction from every individual CATE_model, flag for CATE truncation, average CATE predition used for treatment decision, treatment decision}
#' \item{\code{k_non_na}}{equal to number of k_folds used for nuisance models. Ignore-- placeholder to keep same format as Results object from internal compute estimates}
#' }
compute_estimates_external <- function(df, 
                                       Y_name, 
                                       A_name, 
                                       W_list, 
                                       Z_list, 
                                       CATE_preds, 
                                       threshold, 
                                       nuisance_models,
                                       sl.library.outcome,
                                       sl.library.treatment,
                                       sl.library.missingness,
                                       k_folds,
                                       ps_trunc_level,
                                       outcome_type, 
                                       id_name = NULL){
  
  # Fit nuisance or other model for pseudo-outcome if not pre-fit
  if(is.null(nuisance_models)){
    
    nuisance_output <- drotr::learn_nuisance(df, 
                                             Y_name, 
                                             A_name, 
                                             W_list, 
                                             id_name, 
                                             sl.library.outcome, 
                                             sl.library.treatment,
                                             sl.library.missingness, 
                                             outcome_type, 
                                             k_folds = k_folds, 
                                             ps_trunc_level)
    
    nuisance_models <- nuisance_output$nuisance_models
    
  }
  
  # List to hold results at each threshold
  results_list_threshold <- vector("list", length = length(threshold))
  
  for(t_idx in 1:length(threshold)){
    
    # Extract sign of threshold indicating positive or negative outcome, make threshold numeric
    t <- threshold[t_idx]
    
    sign <- substr(t, 1, 1)
    if(sign != "+" & sign != "-"){ # assume if user did not specify sign, it is positive
      sign <- "+"
    }
    
    t <- as.numeric(t)
    
    # Get treatment decision for avg_CATE
    if(sign == "-"){
      d_pred <- ifelse(CATE_preds$avg_CATE < threshold, 1, 0)
    } else {
      d_pred <- ifelse(CATE_preds$avg_CATE > threshold, 1, 0)
    }
  
    # Get treatment decisions & effect estimates for given nuisance model
    
    k_fold_EY_Ad_dZ1 <- vector("list", length = k_folds )
    k_fold_EY_A0_dZ1 <- vector("list", length = k_folds )
    k_fold_E_dZ1 <- vector("list", length = k_folds )
    k_fold_subgroup_effect <- vector("list", length = k_folds )
    k_fold_treatment_effect <- vector("list", length = k_folds )
    k_fold_inf_fn_matrix <- vector("list", length = k_folds )
    k_fold_subgroup_effect_dZ0 <- vector("list", length = k_folds )
    k_fold_compare_subgroup_effect <- vector("list", length = k_folds )
    
    for(k in 1:k_folds){
      
      nuisance_models_k <- nuisance_models[[k]]
      
      # Estimate treatment effects 
      compute_est_output <- drotr::aiptw_tes(df, 
                                           d_pred, 
                                           Y_name, 
                                           A_name, 
                                           W_list, 
                                           Z_list,
                                           nuisance_models_k,
                                           ps_trunc_level)
      
      k_fold_EY_Ad_dZ1[[k]] <- compute_est_output$EY_Ad_dZ1
      k_fold_EY_A0_dZ1[[k]] <- compute_est_output$EY_A0_dZ1
      k_fold_E_dZ1[[k]] <- compute_est_output$E_dZ1
      k_fold_subgroup_effect[[k]] <- compute_est_output$subgroup_effect
      k_fold_treatment_effect[[k]] <- compute_est_output$treatment_effect
      k_fold_subgroup_effect_dZ0[[k]] <- compute_est_output$subgroup_effect_dZ0
      k_fold_compare_subgroup_effect[[k]] <- compute_est_output$compare_subgroup_effect
      k_fold_inf_fn_matrix[[k]] <- compute_est_output$inf_fn_matrix
      
    }
    
    # Aggregate list of results across folds into dataframe
    k_fold_EY_Ad_dZ1 <- do.call(rbind, k_fold_EY_Ad_dZ1)
    k_fold_EY_A0_dZ1 <- do.call(rbind, k_fold_EY_A0_dZ1)
    k_fold_E_dZ1 <- do.call(rbind, k_fold_E_dZ1)
    k_fold_subgroup_effect <- do.call(rbind, k_fold_subgroup_effect)
    k_fold_subgroup_effect_dZ0 <- do.call(rbind, k_fold_subgroup_effect_dZ0)
    k_fold_treatment_effect <- do.call(rbind, k_fold_treatment_effect)
    k_fold_compare_subgroup_effect <- do.call(rbind, k_fold_compare_subgroup_effect)
    
    # note removed NA folds logic- should not apply when using avg CATE/external dataset
    
    num_obs <- length(d_pred) 
    
    aggregated_results <- data.frame(
      threshold = t,
      aiptw_EY_Ad_dZ1 = mean(k_fold_EY_Ad_dZ1$aiptw),
      se_aiptw_EY_Ad_dZ1 = sqrt((sum(k_fold_EY_Ad_dZ1$var_aug) / k_folds) / num_obs),
      plug_in_est_EY_Ad_dZ1 = mean(k_fold_EY_Ad_dZ1$plug_in_est),
      mean_aug_EY_Ad_dZ1 = mean(k_fold_EY_Ad_dZ1$mean_aug),
      aiptw_EY_A0_dZ1 = mean(k_fold_EY_A0_dZ1$aiptw),
      se_aiptw_EY_A0_dZ1 = sqrt((sum(k_fold_EY_A0_dZ1$var_aug) / k_folds) / num_obs),
      plug_in_est_EY_A0_dZ1 = mean(k_fold_EY_A0_dZ1$plug_in_est),
      mean_aug_EY_A0_dZ1 = mean(k_fold_EY_A0_dZ1$mean_aug),
      E_dZ1 = mean(k_fold_E_dZ1$E_dZ1),
      se_E_dZ1 = sqrt((sum(k_fold_E_dZ1$var_E_dZ1) / k_folds) / num_obs),
      subgroup_effect = mean(k_fold_subgroup_effect$subgroup_effect),
      se_subgroup_effect = sqrt((sum(k_fold_subgroup_effect$var_subgroup_effect) / k_folds) / num_obs),
      subgroup_effect_dZ0 = mean(k_fold_subgroup_effect_dZ0$subgroup_effect_dZ0),
      se_subgroup_effect_dZ0 = sqrt((sum(k_fold_subgroup_effect_dZ0$var_subgroup_effect_dZ0) / k_folds) / num_obs),
      treatment_effect = mean(k_fold_treatment_effect$treatment_effect),
      se_treatment_effect = sqrt((sum(k_fold_treatment_effect$var_treatment_effect) / k_folds) / num_obs),
      compare_subgroup_effect = mean(k_fold_compare_subgroup_effect$compare_subgroup_effect),
      se_compare_subgroup_effect = sqrt((sum(k_fold_compare_subgroup_effect$var_compare_subgroup_effect) / k_folds) / num_obs)
    )
    
    k_fold_results <- list(
      k_fold_EY_Ad_dZ1 = k_fold_EY_Ad_dZ1,
      k_fold_EY_A0_dZ1 = k_fold_EY_A0_dZ1,
      k_fold_E_dZ1 = k_fold_E_dZ1,
      k_fold_subgroup_effect = k_fold_subgroup_effect,
      k_fold_subgroup_effect_dZ0 = k_fold_subgroup_effect_dZ0,
      k_fold_treatment_effect = k_fold_treatment_effect,
      k_fold_compare_subgroup_effect = k_fold_compare_subgroup_effect,
      influence_fns = k_fold_inf_fn_matrix
    )
    
    
    decision_df <- cbind(CATE_preds, data.frame(d_pred = d_pred))
    
    threshold_results <- list(
      aggregated_results = aggregated_results,
      k_fold_results = k_fold_results,
      decision_df = decision_df,
      k_non_na = k_folds
    )
    
    class(threshold_results) <- "Results"
    
    results_list_threshold[[t_idx]] <- threshold_results
    
  }
  
  threshold_names <- paste("threshold = ", threshold)
  names(results_list_threshold) <- threshold_names
  return(results_list_threshold)
  
}

#' Function to estimate treatment effects using AIPTW for external data
#' 
#' @param df dataframe containing external dataset to apply rule(s) to; must contain Z_list variables that are same as pre-trained rule(s)
#' @param d_pred treatment decisions from single CATE model and threshold combo
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param W_list character vector containing names of covariates in the dataframe to be used for fitting nuisance models
#' @param Z_list character vector containing names of variables in df used to fit CATE model (variables used in treatment rule; must be same names as used in pre-fit CATE model(s))
#' @param nuisance_models list of objects of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models 
#' @param ps_trunc_level numeric evel below which propensity scores will be truncated (to avoid errors in computing AIPTW)
#' 
#' @export
#' 
#' @returns list of effect estimates, influence function matrix
aiptw_tes <- function(df, 
                      d_pred, 
                      Y_name, 
                      A_name, 
                      W_list, 
                      Z_list,
                      nuisance_models,
                      ps_trunc_level){
  
  outcome_model <- nuisance_models$outcome_model
  treatment_model <- nuisance_models$treatment_model
  missingness_model <- nuisance_models$missingness_model
  
  # Full set of observations
  I_Y <- ifelse(is.na(df[[Y_name]]), 1, 0) #indicator for Y missing
  Y <- df[[Y_name]]
  A <- df[[A_name]]
  W <- df[, W_list, drop = FALSE]
  Z <- df[, Z_list, drop = FALSE]
  
  # idxes of observations that are recommended treatment (d_pred == 1)
  idx_sub <- which(d_pred == 1)
  idx_sub0 <- which(d_pred == 0) # for effect in d(Z) = 0 subgroup
  
  mean_dZ <- mean(d_pred)
  augmentation_mean_dZ <- d_pred - mean_dZ
  
  mean_dZ0 <- mean(1 - d_pred)
  augmentation_mean_dZ0 <- (1 -d_pred) - mean_dZ0
  
  # E[Y(d) | d(Z) = 1]
  aiptw_a_1 <- drotr::calc_aiptw(a = 1,
                          A = A,
                          A_name = A_name,
                          W = W,
                          Y = Y,
                          I_Y = I_Y,
                          d_pred = d_pred,
                          mean_dZ = mean_dZ,
                          outcome_model = outcome_model,
                          treatment_model = treatment_model,
                          missingness_model = missingness_model,
                          ps_trunc_level = ps_trunc_level,
                          idx_sub = idx_sub)
  
  # E[Y(0) | d(Z) = 1 ]
  aiptw_a_0 <- drotr::calc_aiptw(a = 0,
                          A = A,
                          A_name = A_name,
                          W = W,
                          Y = Y,
                          I_Y = I_Y,
                          d_pred = d_pred,
                          mean_dZ = mean_dZ,
                          outcome_model = outcome_model,
                          treatment_model = treatment_model,
                          missingness_model = missingness_model,
                          ps_trunc_level = ps_trunc_level,
                          idx_sub = idx_sub)
  
  # E[Y(d) | d(Z) = 0]
  aiptw_a_1_dZ0 <- drotr::calc_aiptw(a = 1,
                              A = A,
                              A_name = A_name,
                              W = W,
                              Y = Y,
                              I_Y = I_Y,
                              d_pred = 1 - d_pred,
                              mean_dZ = 1 - mean_dZ,
                              outcome_model = outcome_model,
                              treatment_model = treatment_model,
                              missingness_model = missingness_model,
                              ps_trunc_level = ps_trunc_level,
                              idx_sub = idx_sub0)
  
  # E[Y(0) | d(Z) = 0]
  aiptw_a_0_dZ0 <- drotr::calc_aiptw(a = 0,
                              A = A,
                              A_name = A_name,
                              W = W,
                              Y = Y,
                              I_Y = I_Y,
                              d_pred = 1 - d_pred,
                              mean_dZ = 1 - mean_dZ,
                              outcome_model = outcome_model,
                              treatment_model = treatment_model,
                              missingness_model = missingness_model,
                              ps_trunc_level = ps_trunc_level,
                              idx_sub = idx_sub0)
  
  # E[Y(d) - Y(0)] = E[Y(d) | d(Z) = 1]*P(d(Z) = 1) - E[Y(0) | d(Z) = 1 ]*P(d(Z) = 1)
  treatment_effect <- (aiptw_a_1[['aiptw']] - aiptw_a_0[['aiptw']]) * mean_dZ
  
  # E[Y(1) - Y(0) | d(Z) = 1] = E[Y(1) | d(Z) = 1] - E[Y(0) | d(Z) = 1 ]
  subgroup_effect <- (aiptw_a_1[['aiptw']] - aiptw_a_0[['aiptw']])
  
  # E[Y(1) - Y(0) | d(Z) = 0] = E[Y(1) | d(Z) = 0] - E[Y(0) | d(Z) = 0 ]
  subgroup_effect_dZ0 <- (aiptw_a_1_dZ0[['aiptw']] - aiptw_a_0_dZ0[['aiptw']])
  
  # E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]
  compare_subgroup_effect <- subgroup_effect - subgroup_effect_dZ0
  
  # "n" x 6 matrix
  inf_fn_matrix <- cbind(
    aiptw_a_1$augmentation,
    aiptw_a_0$augmentation,
    augmentation_mean_dZ,
    aiptw_a_1_dZ0$augmentation,
    aiptw_a_0_dZ0$augmentation,
    augmentation_mean_dZ0
  )
  
  # 6x6 covariance matrix
  cov_matrix <- stats::cov(inf_fn_matrix)
  #var_mean_dZ <- as.numeric(stats::var(augmentation_mean_dZ)) should be the same as cov_matrix[3,3]
  
  # Effect in optimally treated subgroup
  gradient_g_subgroup <- matrix(c(
    1, -1, 0, 0, 0, 0
  ), ncol = 1)
  
  var_subgroup_effect <- t(gradient_g_subgroup) %*% cov_matrix %*% gradient_g_subgroup
  
  # Effect in NOT optimally treated subgroup
  gradient_g_subgroup_dZ0 <- matrix(c(
    0, 0, 0, 1, -1, 0
  ), ncol = 1)
  
  var_subgroup_effect_dZ0 <- t(gradient_g_subgroup_dZ0) %*% cov_matrix %*% gradient_g_subgroup_dZ0
  
  # Comparison of optimally and NOT optimally treated subgroups
  gradient_compare_subgroup <- matrix(c(
    1, -1, 0, -1, 1, 0
  ), ncol = 1)
  
  var_compare_subgroup <- t(gradient_compare_subgroup) %*% cov_matrix %*% gradient_compare_subgroup
  
  # Effect overall
  gradient_g <- matrix(c(
    mean_dZ, -mean_dZ, aiptw_a_1[['aiptw']] - aiptw_a_0[['aiptw']], 0, 0, 0
  ) , ncol = 1)
  
  var_treatment_effect <- t(gradient_g) %*% cov_matrix %*% gradient_g
  
  # add columns for the influence functions of subgroup_effect and treatment_effect
  inf_fn_subgroup_effect <- as.numeric(inf_fn_matrix %*% gradient_g_subgroup)
  inf_fn_subgroup_effect_dZ0 <- as.numeric(inf_fn_matrix %*% gradient_g_subgroup_dZ0)
  inf_fn_compare_subgroup <- as.numeric(inf_fn_matrix %*% gradient_compare_subgroup)
  inf_fn_treatment_effect <- as.numeric(inf_fn_matrix %*% gradient_g)
  inf_fn_matrix <- cbind(inf_fn_matrix, inf_fn_subgroup_effect, inf_fn_subgroup_effect_dZ0, inf_fn_treatment_effect)
  
  return(list(EY_Ad_dZ1 = data.frame(aiptw = aiptw_a_1[['aiptw']], plug_in_est = aiptw_a_1[['plug_in_est']],
                                     mean_aug = mean(aiptw_a_1[['augmentation']]), var_aug = aiptw_a_1[['var_aug']]),
              EY_A0_dZ1 = data.frame(aiptw = aiptw_a_0[['aiptw']], plug_in_est = aiptw_a_0[['plug_in_est']],
                                     mean_aug = mean(aiptw_a_0[['augmentation']]), var_aug = aiptw_a_0[['var_aug']]),
              E_dZ1 = data.frame(E_dZ1 = mean_dZ, var_E_dZ1 = cov_matrix[3,3]),
              subgroup_effect = data.frame(subgroup_effect = subgroup_effect, var_subgroup_effect = var_subgroup_effect),
              treatment_effect = data.frame(treatment_effect = treatment_effect, var_treatment_effect = var_treatment_effect),
              subgroup_effect_dZ0 = data.frame(subgroup_effect_dZ0 = subgroup_effect_dZ0, var_subgroup_effect_dZ0 = var_subgroup_effect_dZ0),
              compare_subgroup_effect = data.frame(compare_subgroup_effect = compare_subgroup_effect, var_compare_subgroup_effect = var_compare_subgroup),
              inf_fn_matrix = inf_fn_matrix))# ,
              #df_decisions = df_decisions))
  
}

