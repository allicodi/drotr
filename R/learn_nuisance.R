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

  # SuperLearner(... , control = list(saveCVFitLibrary = TRUE))
  # fit$validRows for assigments
  # get validrows from outcome, then pass validrows into each

  # for v in 1:V
  #   retrive vth trainings sample specific fits from $cvFitLibrary
  #   get predictions that i need to make pseudo-oucome from these models
  #   predict from each of those models on the data in  validRows[[v]]
  # combine predictions from all models into a single prediction using the SL.weights from the fit
  # build ensemble model ourself

  # pass same validrows to the CATE model

  # Full set of observations
  I_Y <- ifelse(is.na(df[[Y_name]]), 1, 0) #indicator for Y missing
  Y <- df[[Y_name]]
  A_vec <- df[[A_name]]
  A <- df[,A_name, drop=FALSE]
  W <- df[,W_list, drop=FALSE]

  # Complete observations only (no NA)
  df_complete <- df[!is.na(df[[Y_name]]), ] #dataframe removing missing Ys
  Y_complete <- df_complete[[Y_name]]
  A_complete <- df_complete[, A_name, drop=FALSE]
  W_complete <- df_complete[, W_list, drop = FALSE]

  ### 1 - Fit outcome model ###

  if(outcome_type == "gaussian"){
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::gaussian(),
                          cvControl = list(V=10), SL.library = sl.library.outcome, control = list(saveCVFitLibrary = TRUE))
  } else {
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::binomial(),
                          cvControl = list(V=10), SL.library = sl.library.outcome, control = list(saveCVFitLibrary = TRUE))
  }

  validRowsComplete <- muhat$validRows

  ### 2 - Fit treatment model ###

  pihat <- SuperLearner::SuperLearner(Y = A_vec, X = W, family = stats::binomial(),
                        cvControl = list(V=10), SL.library = sl.library.treatment)

  validRowsAll <- pihat$validRows

  ### 3 - Fit missingness model ###

  deltahat <- SuperLearner::SuperLearner(Y = I_Y, X = data.frame(A, W), family = stats::binomial(),
                              cvControl = list(V=10, validRows = validRowsAll), SL.library = sl.library.missingness)

  ### 4 - Create pseudo-outcome  ###

  # a. get pred from outcome model under obs trt and obs cov (outcome model)
  om_obs <- stats::predict(muhat, data.frame(A, W), type = 'response')

  # b. predict from trt = 1 and obs cov (outcome model)
  om_trt1_df <- data.frame(1, W)
  names(om_trt1_df)[1] <- A_name
  om_trt1 <- stats::predict(muhat, om_trt1_df, type = 'response')

  # c. predict from trt = 0 and obs cov (outcome model)
  om_trt0_df <- data.frame(0, W)
  names(om_trt0_df)[1] <- A_name
  om_trt0 <- stats::predict(muhat, om_trt0_df, type = 'response')

  # d. predict treatment from obs cov & truncate observations that are too small (treatment model)
  tm_obs_complete <- stats::predict(pihat, data.frame(W), type = 'response')
  tm_obs_complete$pred[tm_obs_complete$pred < ps_trunc_level] <- ps_trunc_level

  # e. predict missingness under obs cov, obs trt (missingness model)
  missingness_obs_complete <- stats::predict(deltahat, data.frame(A, W), type = 'response')
  missingness_obs_complete$pred[missingness_obs_complete$pred < ps_trunc_level] <- ps_trunc_level

  # f. Estimate CATE (pseudo-outcome)
  p1 <- ((2*A_vec - 1)*as.numeric(!I_Y) ) / (tm_obs_complete$pred * (1-missingness_obs_complete$pred))
  p2 <- Y - om_obs$pred
  p2 <- ifelse(is.na(p2), 0, p2)
  p3 <- om_trt1$pred - om_trt0$pred

  CATE_hat <- p1*p2 + p3

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
