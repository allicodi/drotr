#' Make treatment decisions and compute estimated outcomes/treatment effects
#'
#' @param df dataframe containing dataset
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param W_list character vector containing names of covariates in the dataframe to be used in nuisance models
#' @param Z_list character vector containing names of variables in df used to fit CATE model (variables used in treatment rule)
#' @param k_fold_assign_and_CATE dataframe containing ids, fold assignments, and CATE estimate for each observation in df
#' @param CATE_models list of discrete SuperLearner models for CATE from each fold
#' @param nuisance_models list of objects of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models
#' @param threshold character vector of decision thresholds for CATE to determine OTR. Values should be positive if `Y_name` is desirable outcome, negative if `Y_name` is undesirable outcome. If threshold is 0, use +0 for desirable, -0 for undesirable.
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#' @param truncate_CATE logical to indicate if large CATE predictions should be truncated at -1 and 1 (default = TRUE)
#'
#' @import stats
#' @export
#'
#' @returns
#' \describe{
#'  List of objects of class `Results`. Each object contains the following for a given threshold:
#'  \item{\code{aggregated_results}}{dataframe of overall results aggregated across `k` folds for given threshold}
#'  \item{\code{k_fold_results}}{list of results by fold for given threshold}
#'  \item{\code{decision_df}}{original dataset with decision made for each observation at a given threshold}
#'  \item{\code{k_non_na}}{folds that did not have proportion treated = 1 or = 0 (causing some effect estimates to be NA)}
#'  }
compute_estimates <- function(df, Y_name, A_name, W_list, Z_list,
                              k_fold_assign_and_CATE,
                              nuisance_models, CATE_models,
                              threshold, ps_trunc_level = 0.01, truncate_CATE){

  k_folds <- max(k_fold_assign_and_CATE$k)

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

    k_fold_EY_Ad_dZ1 <- vector("list", length = k_folds)
    k_fold_EY_A0_dZ1 <- vector("list", length = k_folds)
    k_fold_E_dZ1 <- vector("list", length = k_folds)
    k_fold_subgroup_effect <- vector("list", length = k_folds)
    k_fold_treatment_effect <- vector("list", length = k_folds)
    k_fold_inf_fn_matrix <- vector("list", length = k_folds)
    k_fold_subgroup_effect_dZ0 <- vector("list", length = k_folds)
    k_fold_compare_subgroup_effect <- vector("list", length = k_folds)
    k_fold_decisions <- data.frame()

    for(k in 1:k_folds){
      # get patient ids that are NOT in the kth fold
      kth_subset_ids <- k_fold_assign_and_CATE$id[k_fold_assign_and_CATE$k == k]

      # estimate CATE hat model kth fold (testing data)
      df_est <- df[!(df$id %in% kth_subset_ids), , drop = FALSE]

      #get models used in kth fold
      nuisance_model <- nuisance_models[[k]]
      CATE_model <- CATE_models[[k]]

      # Find:
      # (1) dataframe with estimated treatment effect among optimally treated
      # (2) dataframe with estimated outcome if not treated among those who should be treated by decision rule
      # (3) overall treatment effect among optimally treated
      # (4) influence function matrix
      # (5) original kth fold data with corresponding treatment decisions
      compute_est_output <- compute_estimate_k(df_est, Y_name, A_name, W_list, Z_list,
                                               CATE_model, nuisance_model,
                                               sign, t, ps_trunc_level, truncate_CATE)

      k_fold_EY_Ad_dZ1[[k]] <- compute_est_output$EY_Ad_dZ1
      k_fold_EY_A0_dZ1[[k]] <- compute_est_output$EY_A0_dZ1
      k_fold_E_dZ1[[k]] <- compute_est_output$E_dZ1
      k_fold_subgroup_effect[[k]] <- compute_est_output$subgroup_effect
      k_fold_treatment_effect[[k]] <- compute_est_output$treatment_effect
      k_fold_subgroup_effect_dZ0[[k]] <- compute_est_output$subgroup_effect_dZ0
      k_fold_compare_subgroup_effect[[k]] <- compute_est_output$compare_subgroup_effect
      k_fold_inf_fn_matrix[[k]] <- compute_est_output$inf_fn_matrix

      kth_decision_df <- compute_est_output$df_decisions

      k_fold_decisions <- rbind(k_fold_decisions, data.frame(id = kth_decision_df$id,
                                                             k = k,
                                                             threshold = t,
                                                             CATE_pred = kth_decision_df$CATE_pred,
                                                             decision = kth_decision_df$d_pred,
                                                             truncated_CATE = kth_decision_df$trunc_flag))
    }

    # Aggregate list of results across folds into dataframe
    k_fold_EY_Ad_dZ1 <- do.call(rbind, k_fold_EY_Ad_dZ1)
    k_fold_EY_A0_dZ1 <- do.call(rbind, k_fold_EY_A0_dZ1)
    k_fold_E_dZ1 <- do.call(rbind, k_fold_E_dZ1)
    k_fold_subgroup_effect <- do.call(rbind, k_fold_subgroup_effect)
    k_fold_subgroup_effect_dZ0 <- do.call(rbind, k_fold_subgroup_effect_dZ0)
    k_fold_treatment_effect <- do.call(rbind, k_fold_treatment_effect)
    k_fold_compare_subgroup_effect <- do.call(rbind, k_fold_compare_subgroup_effect)

    # If any folds contained NA, do not count in computing overall results (cases when a fold recommends treatment to everybody or nobody)
    k_non_na <- which(!(is.na(k_fold_subgroup_effect$var_subgroup_effect) | is.na(k_fold_EY_Ad_dZ1$var_aug)))
  
    non_na_k_fold_EY_Ad_dZ1 <- k_fold_EY_Ad_dZ1[k_non_na,]
    non_na_k_fold_EY_A0_dZ1 <- k_fold_EY_A0_dZ1[k_non_na,]
    non_na_k_fold_E_dZ1 <- k_fold_E_dZ1[k_non_na,]
    non_na_k_fold_subgroup_effect <- k_fold_subgroup_effect[k_non_na,]
    non_na_k_fold_subgroup_effect_dZ0 <- k_fold_subgroup_effect_dZ0[k_non_na,]
    non_na_k_fold_treatment_effect <- k_fold_treatment_effect[k_non_na,]
    non_na_k_fold_compare_subgroup_effect <- k_fold_compare_subgroup_effect[k_non_na,]

    # get number of non-NA folds + observations in them
    k_folds_non_na <- length(k_non_na)
    num_obs <- nrow(k_fold_decisions[k_fold_decisions$k %in% k_non_na,])

    aggregated_results <- data.frame(
      threshold = t,
      aiptw_EY_Ad_dZ1 = mean(non_na_k_fold_EY_Ad_dZ1$aiptw),
      se_aiptw_EY_Ad_dZ1 = sqrt((sum(non_na_k_fold_EY_Ad_dZ1$var_aug) / k_folds_non_na) / num_obs),
      plug_in_est_EY_Ad_dZ1 = mean(non_na_k_fold_EY_Ad_dZ1$plug_in_est),
      mean_aug_EY_Ad_dZ1 = mean(non_na_k_fold_EY_Ad_dZ1$mean_aug),
      aiptw_EY_A0_dZ1 = mean(non_na_k_fold_EY_A0_dZ1$aiptw),
      se_aiptw_EY_A0_dZ1 = sqrt((sum(non_na_k_fold_EY_A0_dZ1$var_aug) / k_folds_non_na) / num_obs),
      plug_in_est_EY_A0_dZ1 = mean(non_na_k_fold_EY_A0_dZ1$plug_in_est),
      mean_aug_EY_A0_dZ1 = mean(non_na_k_fold_EY_A0_dZ1$mean_aug),
      E_dZ1 = mean(non_na_k_fold_E_dZ1$E_dZ1),
      se_E_dZ1 = sqrt((sum(non_na_k_fold_E_dZ1$var_E_dZ1) / k_folds_non_na) / num_obs),
      subgroup_effect = mean(non_na_k_fold_subgroup_effect$subgroup_effect),
      se_subgroup_effect = sqrt((sum(non_na_k_fold_subgroup_effect$var_subgroup_effect) / k_folds_non_na) / num_obs),
      subgroup_effect_dZ0 = mean(non_na_k_fold_subgroup_effect_dZ0$subgroup_effect_dZ0),
      se_subgroup_effect_dZ0 = sqrt((sum(non_na_k_fold_subgroup_effect_dZ0$var_subgroup_effect_dZ0) / k_folds_non_na) / num_obs),
      treatment_effect = mean(non_na_k_fold_treatment_effect$treatment_effect),
      se_treatment_effect = sqrt((sum(non_na_k_fold_treatment_effect$var_treatment_effect) / k_folds_non_na) / num_obs),
      compare_subgroup_effect = mean(non_na_k_fold_compare_subgroup_effect$compare_subgroup_effect),
      se_compare_subgroup_effect = sqrt((sum(non_na_k_fold_compare_subgroup_effect$var_compare_subgroup_effect) / k_folds_non_na) / num_obs)
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

    decision_df <- k_fold_decisions

    threshold_results <- list(
      aggregated_results = aggregated_results,
      k_fold_results = k_fold_results,
      decision_df = decision_df,
      k_non_na = k_non_na
    )

    class(threshold_results) <- "Results"

    results_list_threshold[[t_idx]] <- threshold_results

  }

  threshold_names <- paste("threshold = ", threshold)
  names(results_list_threshold) <- threshold_names
  return(results_list_threshold)

}

#' Make treatment decisions and compute estimated outcomes/treatment effects in kth fold
#'
#' @param df dataframe containing testing dataset
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param W_list character vector containing names of covariates in the dataframe to be used in nuisance models
#' @param Z_list character vector containing names of variables in df used to fit CATE model (variables used in treatment rule)
#' @param CATE_model discrete SuperLearner model for CATE
#' @param nuisance object of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models
#' @param sign sign of threshold, + indicates `Y` is a desirable outcome, - indicates `Y` is an undesirable outcome
#' @param threshold character vector of decision thresholds for CATE to determine OTR. Values should be positive if `Y_name` is desirable outcome, negative if `Y_name` is undesirable outcome. If threshold is 0, use +0 for desirable, -0 for undesirable.
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#' @param truncate_CATE logical to indicate if large CATE predictions should be truncated at -1 and 1 (default = TRUE)
#'
#' @returns
#' \describe{
#'        (1) dataframe with E[Y(1) | d(Z) = 1] -- estimated treatment effect among optimally treated
#'        (2) dataframe with E[Y(0) | d(Z) = 1] -- estimated outcome if not treated among those who should be treated by decision rule
#'        (3) dataframe with E[d(Z) = 1] -- estimated proportion treated under optimal treatment rule
#'        (4) dataframe with E[Y(1) - Y(0) | d(Z) = 1] -- estimated subgroup effect
#'        (5) dataframe with E[Y(d) - Y(0) ] -- overall treatment effect among optimally treated
#'        (6) dataframe with E[Y(1) - Y(0) | d(Z) = 0] -- estimated subgroup effect in the untreated
#'        (7) dataframe with E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0] -- comparison of subgroup effects between those recommended and not recommended treatment
#'        (8) influence function matrix
#'        (9) original kth fold data with corresponding treatment decisions
#'  }
#'
#' @keywords internal
compute_estimate_k <- function(df, Y_name, A_name, W_list, Z_list,
                             CATE_model, nuisance,
                             sign, threshold, ps_trunc_level = 0.01, truncate_CATE){

  # extract models from nuisance
  outcome_model <- nuisance$outcome_model
  treatment_model <- nuisance$treatment_model
  missingness_model <- nuisance$missingness_model

  # Full set of observations
  I_Y <- ifelse(is.na(df[[Y_name]]), 1, 0) #indicator for Y missing
  Y <- df[[Y_name]]
  A <- df[[A_name]]
  W <- df[, W_list, drop = FALSE]
  Z <- df[, Z_list, drop = FALSE]

  ### Step 1: using df_est, get prediction from CATE model and find observations that meet treatment threshold
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

  if(sign == "-"){
    d_pred <- ifelse(CATE_pred < threshold, 1, 0)
  } else {
    d_pred <- ifelse(CATE_pred > threshold, 1, 0)
  }

  # idxes of observations that are recommended treatment (d_pred == 1)
  idx_sub <- which(d_pred == 1)
  idx_sub0 <- which(d_pred == 0) # for effect in d(Z) = 0 subgroup

  # add decisions to kth fold dataframe to return later on
  df_decisions <- cbind(df, d_pred, data.frame(CATE_pred = CATE_pred,
                                               trunc_flag = trunc_flag))

  ### Step 2: Find P(d(Z) = 1) by taking mean(d(Z)) created in step 1
  mean_dZ <- mean(d_pred)

  # E[Y(d) | d(Z) = 1]
  aiptw_a_1 <- calc_aiptw(a = 1,
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
  aiptw_a_0 <- calc_aiptw(a = 0,
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
  aiptw_a_1_dZ0 <- calc_aiptw(a = 1,
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
  aiptw_a_0_dZ0 <- calc_aiptw(a = 0,
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

  ### Step 10: estimate of influence function & its variance
  mean_dZ <- mean(d_pred)
  augmentation_mean_dZ <- d_pred - mean_dZ

  mean_dZ0 <- mean(1 - d_pred)
  augmentation_mean_dZ0 <- (1 -d_pred) - mean_dZ0

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

  # return  (1) dataframe with E[Y(1) | d(Z) = 1] -- estimated treatment effect among optimally treated
  #         (2) dataframe with E[Y(0) | d(Z) = 1] -- estimated outcome if not treated among those who should be treated by decision rule
  #         (3) dataframe with E[d(Z) = 1] -- estimated proportion treated under optimal treatment rule
  #         (4) dataframe with E[Y(1) - Y(0) | d(Z) = 1] -- estimated subgroup effect
  #         (5) dataframe with E[Y(d) - Y(0) ] -- overall treatment effect among optimally treated
  #         (6) dataframe with E[Y(1) - Y(0) | d(Z) = 0] -- estimated subgroup effect in the untreated
  #         (7) dataframe with E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0] -- comparison of subgroup effects between those recommended and not recommended treatment
  #         (8) influence function matrix
  #         (9) original kth fold data with corresponding treatment decisions
  return(list(EY_Ad_dZ1 = data.frame(aiptw = aiptw_a_1[['aiptw']], plug_in_est = aiptw_a_1[['plug_in_est']],
                         mean_aug = mean(aiptw_a_1[['augmentation']]), var_aug = aiptw_a_1[['var_aug']]),
              EY_A0_dZ1 = data.frame(aiptw = aiptw_a_0[['aiptw']], plug_in_est = aiptw_a_0[['plug_in_est']],
                         mean_aug = mean(aiptw_a_0[['augmentation']]), var_aug = aiptw_a_0[['var_aug']]),
              E_dZ1 = data.frame(E_dZ1 = mean_dZ, var_E_dZ1 = cov_matrix[3,3]),
              subgroup_effect = data.frame(subgroup_effect = subgroup_effect, var_subgroup_effect = var_subgroup_effect),
              treatment_effect = data.frame(treatment_effect = treatment_effect, var_treatment_effect = var_treatment_effect),
              subgroup_effect_dZ0 = data.frame(subgroup_effect_dZ0 = subgroup_effect_dZ0, var_subgroup_effect_dZ0 = var_subgroup_effect_dZ0),
              compare_subgroup_effect = data.frame(compare_subgroup_effect = compare_subgroup_effect, var_compare_subgroup_effect = var_compare_subgroup),
              inf_fn_matrix = inf_fn_matrix,
              df_decisions = df_decisions))
}

#' Function to calculate AIPTW
#'
#' @param a value of treatment to set for all observations
#' @param A vector of observed treatments
#' @param A_name name of treatment variable
#' @param W dataframe of covariates
#' @param Y vector of observed outcome variable
#' @param I_Y indicator that Y is missing (I_Y = 1 if missing, I_Y = 0 if observed)
#' @param d_pred vector of treatment decisions (d_pred = 1 if recommend treatment under rule)
#' @param mean_dZ proportion recommended treatment under rule
#' @param outcome_model outcome Nuisance model
#' @param treatment_model treatment Nuisance model
#' @param missingness_model missingness Nuisance model
#' @param ps_trunc_level threshold to truncate propensity score
#' @param idx_sub indices of observations that are recommended treatment (d_pred == 1)
#'
#' @export
#'
#' @returns
#' \describe{
#'         (1) a - value of treatment
#'         (2) aiptw - vector aiptw estimators
#'         (3) plug_in_est - vector of plug-in estimators
#'         (4) augmentation - vector of augmentation terms
#'         (5) var_aug - variance of augmentation term
#'         }
calc_aiptw <- function(a, A, A_name, W, Y, I_Y,
                       d_pred, mean_dZ, outcome_model, treatment_model, missingness_model,
                       ps_trunc_level, idx_sub){

  ### Step 3: get 1-prediction from missingness model for everyone in df_est P(Delta = 1 | ...)
  miss_pred_output_df <- data.frame(a, W)
  names(miss_pred_output_df)[1] <- A_name

  miss_pred_output <- stats::predict(missingness_model, miss_pred_output_df, type = 'response')
  miss_pred <- 1 - miss_pred_output$pred # proportion of non-missing observations

  ### Step 4: estimate P(A = a | W)
  tm_given_rec <- stats::predict(treatment_model, W, type = 'response')
  tm_given_rec$pred <- ifelse(tm_given_rec$pred < ps_trunc_level, ps_trunc_level, tm_given_rec$pred) #truncate if too small

  if(a == 0){
    tm_given_rec$pred <- 1 - tm_given_rec$pred
  }

  ### Step 5: predict from outcome model (muhat) on df_est, setting A = a
  om_trta_df <- data.frame(a, W)
  names(om_trta_df)[1] <- A_name
  om_trta <- stats::predict(outcome_model, om_trta_df, type = 'response')

  ### Step 6: "compute plug-in estimate"
  plug_in_est <- mean(om_trta$pred[idx_sub])

  ### Step 7: compute augmentation term
  I_A1_notMiss_d1 <- ifelse(A == a & I_Y == 0 & d_pred == 1, 1, 0)
  denom <- miss_pred * tm_given_rec$pred * mean_dZ
  Y_minus_om_trta <- Y - om_trta$pred
  Y_minus_om_trta <- ifelse(is.na(Y_minus_om_trta), 0, Y_minus_om_trta) # will get 0'd out by indicator, avoid error
  ind_over_prob <- d_pred / mean_dZ
  om_trta_minus_plugin <- om_trta$pred - plug_in_est

  augmentation <- (I_A1_notMiss_d1 / denom * Y_minus_om_trta) + (ind_over_prob*om_trta_minus_plugin)

  ### Step 8: compute aiptw estimate = plug-in estimate from 6 + mean(augmentation term from step 7)
  aiptw <- plug_in_est + mean(augmentation)

  ### Step 9: compute variance of augmentation (used for standard error at end)
  var_aug <- as.numeric(stats::var(augmentation))

  return(list(
    a = a,
    aiptw = aiptw,
    plug_in_est = plug_in_est,
    augmentation = augmentation,
    var_aug = var_aug
  ))
}

