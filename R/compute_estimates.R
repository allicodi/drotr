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
#' @param threshold decision threshold for CATE to determine OTR. `treatment` should be positive if `Y_name` is desirable outcome, negative if `Y_name` is undesirable outcome
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#'
#' @importFrom dplyr left_join
#' @import stats
#' @export
#'
#' @returns
#' \describe{
#'  \item{\code{overall_results}}{dataframe of overall results aggregated across `k` folds}
#'  \item{\code{EY_A1_d1}}{dataframe of AIPTW for optimally treated in each fold}
#'  \item{\code{EY_A0_d1}}{dataframe of AIPTW for not treating those who should be treated under decision rule in each fold}
#'  \item{\code{treatment_effect}}{dataframe of treatment effect in each fold}
#'  \item{\code{influence_fns}}{influence functions used in each fold}
#'  \item{\code{decision_df}}{original dataset with decision made for each observation}
#'  }
compute_estimates <- function(df, Y_name, A_name, W_list, Z_list,
                              k_fold_assign_and_CATE,
                              nuisance_models, CATE_models,
                              threshold, ps_trunc_level = 0.01){

  # if pid not in df, make pid index
  # if (!"pid" %in% colnames(df)) {
  #   df$pid <- as.numeric(rownames(df))
  # }

  k_folds <- max(k_fold_assign_and_CATE$k)

  k_fold_output_aiptw_a_1 <- list()
  k_fold_output_aiptw_a_0 <- list()
  k_fold_treatment_effect <- list()
  k_fold_inf_fn_matrix <- list()
  k_fold_decisions <- data.frame()

  for(k in 1:k_folds){
    # get patient ids that are in the kth fold
    kth_subset_ids <- k_fold_assign_and_CATE$id[k_fold_assign_and_CATE$k == k]

    # estimate CATE hat model kth fold (testing data)
    df_est <- df[df$id %in% kth_subset_ids, , drop = FALSE]

    # join CATEhat with data by id
    df_est <- dplyr::left_join(df_est, k_fold_assign_and_CATE, by = "id")

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
                                     threshold, ps_trunc_level)

    k_fold_output_aiptw_a_1 <- c(k_fold_output_aiptw_a_1, list(compute_est_output[[1]]))
    k_fold_output_aiptw_a_0 <- c(k_fold_output_aiptw_a_0, list(compute_est_output[[2]]))
    k_fold_treatment_effect <- c(k_fold_treatment_effect, list(compute_est_output[[3]]))
    k_fold_inf_fn_matrix <- c(k_fold_inf_fn_matrix, list(compute_est_output[[4]]))

    kth_decision_df <- compute_est_output[[5]]
    k_fold_decisions <- rbind(k_fold_decisions, data.frame(id = kth_decision_df$id,
                                                           k = k,
                                                           decision = kth_decision_df$d_pred))
  }

  # Aggregate list of results across folds into dataframe
  k_fold_output_aiptw_a_1 <- do.call(rbind, k_fold_output_aiptw_a_1)
  k_fold_output_aiptw_a_0 <- do.call(rbind, k_fold_output_aiptw_a_0)
  k_fold_treatment_effect <- do.call(rbind, k_fold_treatment_effect)

  # If any folds were NA, do not count in computing overall results
  # k_non_na will be the same for aiptw_a_1 and aiptw_a_0, find once
  k_non_na <- which(!is.na(k_fold_output_aiptw_a_1$var_aug))
  k_folds_non_na <- length(k_non_na)
  num_obs <- nrow(k_fold_assign_and_CATE[k_fold_assign_and_CATE$k %in% k_non_na,])

  results <- data.frame(aiptw_a_1 = mean(k_fold_output_aiptw_a_1$aiptw, na.rm=TRUE),
                        se_aiptw_a_1 = sqrt((sum(k_fold_output_aiptw_a_1$var_aug, na.rm=TRUE) / k_folds_non_na) / num_obs),
                        plug_in_est_1 = mean(k_fold_output_aiptw_a_1$plug_in_est, na.rm=TRUE),
                        mean_aug_1 = mean(k_fold_output_aiptw_a_1$mean_aug, na.rm=TRUE),
                        aiptw_a_0 = mean(k_fold_output_aiptw_a_0$aiptw, na.rm=TRUE),
                        se_aiptw_a_0 = sqrt((sum(k_fold_output_aiptw_a_0$var_aug, na.rm=TRUE) / k_folds_non_na) / num_obs),
                        plug_in_est_0 = mean(k_fold_output_aiptw_a_0$plug_in_est, na.rm=TRUE),
                        mean_aug_0 = mean(k_fold_output_aiptw_a_0$mean_aug, na.rm=TRUE),
                        treatment_effect = mean(k_fold_treatment_effect$treatment_effect, na.rm=TRUE),
                        se_treatment_effect = sqrt((sum(k_fold_treatment_effect$variance_treatment_effect, na.rm = TRUE) / k_folds_non_na) / num_obs))

  #Return overall results,
  # results from each fold for aiptw_a_1,
  # results from each fold for aiptw_a_0,
  # results from each fold for treatment effect,
  # influence functions from each fold
  # dataset with id/predictions
  return(list(overall_results = results,
              EY_A1_d1 = k_fold_output_aiptw_a_1,
              EY_A0_d1 = k_fold_output_aiptw_a_0,
              treatment_effect = k_fold_treatment_effect,
              influence_fns = k_fold_inf_fn_matrix,
              decision_df = k_fold_decisions))
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
#' @param threshold decision threshold for CATE to determine OTR. `treatment` should be positive if `Y_name` is desirable outcome, negative if `Y_name` is undesirable outcome
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#'
#' @returns a list with (1) dataframe with estimated treatment effect among optimally treated, (2) dataframe with estimated outcome if not treated among those who should be treated by decision rule, (3) overall treatment effect among optimally treated, (4) influence function matrix, (5) original kth fold data with corresponding treatment decisions
#'
#' @keywords internal
compute_estimate_k <- function(df, Y_name, A_name, W_list, Z_list,
                             CATE_model, nuisance,
                             threshold, ps_trunc_level = 0.01){

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
  CATE_pred <- stats::predict(CATE_model, Z, type = 'response', family = stats::gaussian())

  # d_pred <- ifelse(CATE_pred < -abs(threshold), 1, 0)
  if(threshold < 0){
    d_pred <- ifelse(CATE_pred < threshold, 1, 0)
  } else {
    d_pred <- ifelse(CATE_pred > threshold, 1, 0)
  }

  # idxes of observations that are recommended treatment (d_pred == 1)
  idx_sub <- which(d_pred == 1)

  # add decisions to kth fold dataframe to return later on
  df_decisions <- cbind(df, d_pred)

  ### Step 2: Find P(d(Z) = 1) by taking mean(d(Z)) created in step 1
  mean_dZ <- mean(d_pred)

  #
  # Helper function to calculate AIPTW
  #
  calc_aiptw <- function(a, mean_dZ, outcome_model, treatment_model, missingness_model){

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

  # Call calc_aiptw function for Y(d) and Y(0)

  # E[Y(d) | d(Z) = 1]
  aiptw_a_1 <- calc_aiptw(1, mean_dZ, outcome_model, treatment_model, missingness_model)

  # E[Y(0) | d(Z) = 1 ]
  aiptw_a_0 <- calc_aiptw(0, mean_dZ, outcome_model, treatment_model, missingness_model)

  ### Step 10: estimate of influence function & its variance
  mean_dZ <- mean(d_pred)
  augmentation_mean_dZ <- d_pred - mean_dZ
  var_mean_dZ <- as.numeric(stats::var(augmentation_mean_dZ))

  # E[Y(d) - Y(0)] = E[Y(d) | d(Z) = 1]*P(d(Z) = 1) - E[Y(0) | d(Z) = 1 ]*P(d(Z) = 1)
  treatment_effect <- (aiptw_a_1[['aiptw']] - aiptw_a_0[['aiptw']]) * mean_dZ

  #todo- return this without multiplying by the proportion

  # "n" x 3 matrix
  inf_fn_matrix <- cbind(
    aiptw_a_1$augmentation,
    aiptw_a_0$augmentation,
    augmentation_mean_dZ
  )

  # 3x3 covariance matrix
  cov_matrix <- stats::cov(inf_fn_matrix)
  gradient_g <- matrix(c(
    mean_dZ, -mean_dZ, aiptw_a_1[['aiptw']] - aiptw_a_0[['aiptw']]
  ) , ncol = 1)

  var_treatment_effect <- t(gradient_g) %*% cov_matrix %*% gradient_g

  # add a column for the influence function of treatment_effect
  inf_fn_treatment_effect <- as.numeric(inf_fn_matrix %*% gradient_g)
  inf_fn_matrix <- cbind(inf_fn_matrix, inf_fn_treatment_effect)

  # return  (1) dataframe with estimated treatment effect among optimally treated
  #         (2) dataframe with estimated outcome if not treated among those who should be treated by decision rule
  #         (3) overall treatment effect among optimally treated
  #         (4) influence function matrix
  #         (5) original kth fold data with corresponding treatment decisions
  return(list(data.frame(aiptw = aiptw_a_1[['aiptw']], plug_in_est = aiptw_a_1[['plug_in_est']],
                         mean_aug = mean(aiptw_a_1[['augmentation']]), var_aug = aiptw_a_1[['var_aug']]),
              data.frame(aiptw = aiptw_a_0[['aiptw']], plug_in_est = aiptw_a_0[['plug_in_est']],
                         mean_aug = mean(aiptw_a_0[['augmentation']]), var_aug = aiptw_a_0[['var_aug']]),
              data.frame(treatment_effect = treatment_effect, variance_treatment_effect = var_treatment_effect),
              inf_fn_matrix,
              df_decisions))

}
