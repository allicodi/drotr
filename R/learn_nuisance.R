#' Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats using k-fold cross validation
#'
#' @param df dataframe containing full dataset
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param W_list character vector containing names of covariates in the dataframe to be used for fitting nuisance models
#' @param id_name name of patient id variable in dataset if applicable, will default to NULL and use observation index
#' @param sl.library.outcome character vector of SuperLearner libraries to use to fit the outcome models
#' @param sl.library.treatment character vector of SuperLearner libraries to use to fit the treatment models
#' @param sl.library.missingness character vector of SuperLearner libraries to use to fit the missingness models
#' @param k_folds number of folds for k_fold cross validation
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#' @param outcome_type specifying continuous (outcome_type = "gaussian") or binary (outcome_type = "binomial") outcome Y
#'
#' @import stats
#' @import SuperLearner
#' @export
#'
#' @return
#' \describe{
#'  \item{\code{k_fold_nuisance}}{list of `Nuisance` objects (fit nuisance models) for each of k folds}
#'  \item{\code{CATE_hat}}{dataframe of CATE estimates and fold assignments for each observation}
#'  }
learn_nuisance <- function(df,
                           Y_name,
                           A_name,
                           W_list,
                           id_name = NULL,
                           sl.library.outcome,
                           sl.library.treatment,
                           sl.library.missingness,
                           outcome_type,
                           k_folds = 2,
                           ps_trunc_level = 0.01){

  # if id_name not specified, make index id
  if(is.null(id_name)){
    df$id <- rownames(df) #fix this, don't use rownames, just assign id
  } else {
    df$id <- df[[id_name]]
  }

  # split into k equal sized folds and randomly shuffle
  folds <- cut(seq(1,nrow(df)),breaks=k_folds,labels=FALSE)
  folds <- sample(folds)

  fold_assignments <- data.frame(id = df$id, fold = folds)

  k_fold_nuisance <- vector(mode = "list", length = k_folds)
  k_fold_assign_and_CATE <- data.frame()

  validRowsCompleteByFold <- vector(mode = "list", length = k_folds)
  validRowsAllByFold <- vector(mode = "list", length = k_folds)

  for(k in 1:k_folds){

    df_learn <- df[-which(folds == k),] # train on all but kth fold

    output_learn <- learn_nuisance_k(df_learn,
                                      Y_name,
                                      A_name,
                                      W_list,
                                      sl.library.outcome,
                                      sl.library.treatment,
                                      sl.library.missingness,
                                      outcome_type,
                                      ps_trunc_level)

    k_fold_nuisance[[k]] <- output_learn$nuisance_models

    # note here k is actually denoting the fold that was left out of the training set
    k_fold_assign_and_CATE <- rbind(k_fold_assign_and_CATE, data.frame(id = as.numeric(df_learn$id),
                                                                       k = k,
                                                                       pseudo_outcome = output_learn$pseudo_outcome))
    validRowsCompleteByFold[[k]] <- output_learn$validRowsComplete
    validRowsAllByFold[[k]] <- output_learn$validRowsAll
  }

  return(list(nuisance_models = k_fold_nuisance,
              k_fold_assign_and_CATE = k_fold_assign_and_CATE,
              validRowsComplete = validRowsCompleteByFold,
              validRowsAll = validRowsAllByFold,
              fold_assignments = fold_assignments))
}

#' Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats for kth fold
#'
#' Note training dataset contains all observations EXCEPT those assigned to kth fold
#'
#' @param df dataframe containing training dataset
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param W_list character vector containing names of covariates in the dataframe to be used for fitting nuisance models
#' @param sl.library.outcome character vector of SuperLearner libraries to use to fit the outcome model
#' @param sl.library.treatment character vector of SuperLearner libraries to use to fit the treatment model
#' @param sl.library.missingness character vector of SuperLearner libraries to use to fit the missingness model
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#' @param outcome_type specifying continuous (gaussian) or binomial outcome Y
#'
#' @return
#' \describe{
#'  \item{\code{k_fold_nuisance}}{object of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models}
#'  \item{\code{k_fold_assign_and_CATE}}{numeric vector of CATE estimates for data in df}
#'  }
#' @keywords internal
learn_nuisance_k <- function(df, Y_name, A_name, W_list,
                         sl.library.outcome,
                         sl.library.treatment,
                         sl.library.missingness,
                         outcome_type,
                         ps_trunc_level = 0.01){

  idx_na_Y <- which(is.na(df[[Y_name]]))
  idx_no_na_Y <- which(!is.na(df[[Y_name]]))
  nNA <- length(idx_na_Y)
  n <- nrow(df)

  df_sort <- df[c(idx_na_Y, idx_no_na_Y), ]

  # Full set of observations
  I_Y <- ifelse(is.na(df_sort[[Y_name]]), 1, 0) #indicator for Y missing
  Y <- df_sort[[Y_name]]
  A_vec <- df_sort[[A_name]]
  A <- df_sort[,A_name, drop=FALSE]
  W <- df_sort[,W_list, drop=FALSE]

  # Complete observations only (no NA)
  df_complete <- df_sort[!is.na(df[[Y_name]]), ] #dataframe removing missing Ys
  Y_complete <- df_complete[[Y_name]]
  A_complete <- df_complete[, A_name, drop=FALSE]
  W_complete <- df_complete[, W_list, drop = FALSE]

  if(nNA > 0){
    rows_na_Y <- split(sample(1:nNA), rep(1:10, length = nNA))
  }else{
    rows_na_Y <- vector(mode = "list", length = 0L)
  }
  rows_no_na_Y <- split(sample((nNA+1):n), rep(1:10, length = n-nNA))
  master_validRows <- vector(mode = "list", length = 10)
  for (vv in seq_len(10)) {
    if(length(rows_na_Y) >= vv & length(rows_no_na_Y) >= vv){
      master_validRows[[vv]] <- c(rows_na_Y[[vv]], rows_no_na_Y[[vv]])
    }else if(length(rows_na_Y) < vv){
      master_validRows[[vv]] <- rows_no_na_Y[[vv]]
    }else if(length(rows_no_na_Y) < vv){
      master_validRows[[vv]] <- rows_na_Y[[vv]]
    }
  }

  ### 1 - Fit outcome model ###
  muhat_validRows <- lapply(master_validRows, function(x){
    x - nNA
  })

  if(outcome_type == "gaussian"){
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::gaussian(),
                          cvControl = list(validRows = muhat_validRows), SL.library = sl.library.outcome, control = list(saveCVFitLibrary = TRUE))
  } else {
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::binomial(),
                          cvControl = list(validRows = muhat_validRows), SL.library = sl.library.outcome, control = list(saveCVFitLibrary = TRUE))
  }

  muhat.cvFitLibrary <- muhat$cvFitLibrary
  muhat.coef <- muhat$coef

  ### 2 - Fit treatment model ###

  pihat <- SuperLearner::SuperLearner(Y = A_vec, X = W, family = stats::binomial(),
                        cvControl = list(validRows = master_validRows), SL.library = sl.library.treatment, control = list(saveCVFitLibrary = TRUE))

  pihat.cvFitLibrary <- pihat$cvFitLibrary
  pihat.coef <- pihat$coef

  ### 3 - Fit missingness model ###

  deltahat <- SuperLearner::SuperLearner(Y = I_Y, X = data.frame(A, W), family = stats::binomial(),
                              cvControl = list(validRows = master_validRows), SL.library = sl.library.missingness, control = list(saveCVFitLibrary = TRUE))

  deltahat.cvFitLibrary <- deltahat$cvFitLibrary
  deltahat.coef <- deltahat$coef

  ### 4 - Create pseudo-outcome  ###

  # for v in 1:V
  #   retrive vth trainings sample specific fits from $cvFitLibrary
  #   get predictions that i need to make pseudo-oucome from these models
  #   predict from each of those models on the data in  validRows[[v]]
  # combine predictions from all models into a single prediction using the SL.weights from the fit
  # build ensemble model ourself

  # pass same validrows to the CATE model

  muhat_obs_matrix <- matrix(NA, nrow = nrow(df), ncol = 10)
  muhat_1_matrix <- matrix(NA, nrow = nrow(df), ncol = 10)
  muhat_0_matrix <- matrix(NA, nrow = nrow(df), ncol = 10)
  pihat_matrix <- matrix(NA, nrow = nrow(df), ncol = 10)
  deltahat_matrix <- matrix(NA, nrow = nrow(df), ncol = 10)

  muhat_individ_mod <- matrix(NA, nrow = nrow(df), ncol = length(muhat.v))

  for(v in 1:10){
    muhat.v <- muhat.cvFitLibrary[[v]]
    pihat.v <- pihat.cvFitLibrary[[v]]
    deltahat.v <- deltahat.cvFitLibrary[[v]]

    # should I still be predicting on all of them? or just the ones in the vth fold
    # also muhat was fit using just complete, but we're predicting using all of them?

    # a. get pred from outcome model under obs trt and obs cov (outcome model)
    for(model in 1:length(muhat.v)){
      muhat.pred <- stats::predict(
        muhat.v[[model]],
        newdata = data.frame(A, W)[master_validRows[[v]], , drop = FALSE]
      )
      muhat_individ_mod[master_validRows[[v]], model] <- muhat.pred
    }
#
#     om_obs.v <- muhat_individ_mod %*% muhat.coef
#     muhat_obs_matrix[,v] <- om_obs.v


    # b. predict from trt = 1 and obs cov (outcome model)
    om_trt1_df <- data.frame(1, W)
    names(om_trt1_df)[1] <- A_name

    muhat_individ_mod1 <- matrix(NA, nrow = nrow(df), ncol = length(muhat.v))
    for(model in 1:length(muhat.v)){
      muhat1.pred <- stats::predict(muhat.v[[model]], om_trt1_df, type = 'response')
      muhat_individ_mod1[,model] <- muhat1.pred
    }

    om_trt1.v <- muhat_individ_mod1 %*% muhat.coef
    muhat_1_matrix[,v] <- om_trt1.v

    # c. predict from trt = 0 and obs cov (outcome model)
    om_trt0_df <- data.frame(0, W)
    names(om_trt0_df)[1] <- A_name

    muhat_individ_mod0 <- matrix(NA, nrow = nrow(df), ncol = length(muhat.v))
    for(model in 1:length(muhat.v)){
      muhat0.pred <- stats::predict(muhat.v[[model]], om_trt0_df, type = 'response')
      muhat_individ_mod0[,model] <- muhat0.pred
    }

    om_trt0.v <- muhat_individ_mod0 %*% muhat.coef
    muhat_0_matrix[,v] <- om_trt0.v

    # d. predict treatment from obs cov & truncate observations that are too small (treatment model)

    pihat_individ_mod <- matrix(NA, nrow = nrow(df), ncol = length(pihat.v))
    for(model in 1:length(pihat.v)){
      pihat.pred <- stats::predict(pihat.v[[model]], data.frame(W), type = 'response')
      pihat_individ_mod[,model] <- pihat.pred
    }

    pihat.v <- pihat_individ_mod %*% pihat.coef
    pihat.v[pihat.v < ps_trunc_level] <- ps_trunc_level
    pihat_matrix[,v] <- pihat.v

    # e. predict missingness under obs cov, obs trt (missingness model)
    delta_individ_mod <- matrix(NA, nrow = nrow(df), ncol=length(deltahat.v))
    for(model in 1:length(deltahat.v)){
      deltahat.pred <- stats::predict(deltahat.v[[model]], data.frame(A, W), type = 'response')
      delta_individ_mod[,model] <- deltahat.pred
    }

    deltahat.v <- delta_individ_mod %*% deltahat.coef
    deltahat.v[deltahat.v <- ps_trunc_level] <- ps_trunc_level
    deltahat_matrix[,v] <- deltahat.v

  }

  avg_muhat_obs <- rowMeans(muhat_obs_matrix)
  avg_muhat_1 <- rowMeans(muhat_1_matrix)
  avg_muhat_0 <- rowMeans(muhat_0_matrix)
  avg_pihat_obs <- rowMeans(pihat_matrix)
  avg_deltahat_obs <- rowMeans(deltahat_matrix)

  # f. Estimate CATE (pseudo-outcome)
  p1 <- ((2*A_vec - 1)*as.numeric(!I_Y) ) / (avg_pihat_obs * (1-avg_deltahat_obs))
  p2 <- Y - avg_muhat_obs
  p2 <- ifelse(is.na(p2), 0, p2)
  p3 <- avg_muhat_1 - avg_muhat_0

  CATE_hat <- p1*p2 + p3

#
#
# # ------------------------------------------------------------------------
#
#   # a. get pred from outcome model under obs trt and obs cov (outcome model)
#   om_obs <- stats::predict(muhat, data.frame(A, W), type = 'response')
#
#   # b. predict from trt = 1 and obs cov (outcome model)
#   om_trt1_df <- data.frame(1, W)
#   names(om_trt1_df)[1] <- A_name
#   om_trt1 <- stats::predict(muhat, om_trt1_df, type = 'response')
#
#   # c. predict from trt = 0 and obs cov (outcome model)
#   om_trt0_df <- data.frame(0, W)
#   names(om_trt0_df)[1] <- A_name
#   om_trt0 <- stats::predict(muhat, om_trt0_df, type = 'response')
#
#   # d. predict treatment from obs cov & truncate observations that are too small (treatment model)
#   tm_obs_complete <- stats::predict(pihat, data.frame(W), type = 'response')
#   tm_obs_complete$pred[tm_obs_complete$pred < ps_trunc_level] <- ps_trunc_level
#
#   # e. predict missingness under obs cov, obs trt (missingness model)
#   missingness_obs_complete <- stats::predict(deltahat, data.frame(A, W), type = 'response')
#   missingness_obs_complete$pred[missingness_obs_complete$pred < ps_trunc_level] <- ps_trunc_level
#
#   # f. Estimate CATE (pseudo-outcome)
#   p1 <- ((2*A_vec - 1)*as.numeric(!I_Y) ) / (tm_obs_complete$pred * (1-missingness_obs_complete$pred))
#   p2 <- Y - om_obs$pred
#   p2 <- ifelse(is.na(p2), 0, p2)
#   p3 <- om_trt1$pred - om_trt0$pred
#
#   CATE_hat <- p1*p2 + p3

  # Create Nuisance object
  learned_models <- list(outcome_model = muhat,
                         treatment_model = pihat,
                         missingness_model = deltahat)

  class(learned_models) <- "Nuisance"

  return(list(nuisance_models = learned_models,
              pseudo_outcome = CATE_hat,
              validRowsComplete = validRowsComplete,
              validRowsAll = validRowsAll))
}
