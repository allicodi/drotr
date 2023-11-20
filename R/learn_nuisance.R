#' Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats using k-fold cross validation
#'
#' @param df dataframe containing full dataset
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param sl.library.outcome character vector of SuperLearner libraries to use to fit the outcome models
#' @param sl.library.treatment character vector of SuperLearner libraries to use to fit the treatment models
#' @param sl.library.missingness character vector of SuperLearner libraries to use to fit the missingness models
#' @param k_folds number of folds for k_fold cross validation
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#' @param outcome_type specifying continuous (outcome_type = "gaussian") or binary (outcome_type = "binomial") outcome Y
#'
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
                           sl.library.outcome,
                           sl.library.treatment,
                           sl.library.missingness,
                           outcome_type,
                           k_folds = 2,
                           ps_trunc_level = 0.01){

  # split into k equal sized folds and randomly shuffle
  folds <- cut(seq(1,nrow(df)),breaks=k_folds,labels=FALSE)
  folds <- sample(folds)

  k_fold_nuisance <- vector(mode = "list", length = k_folds)
  k_fold_assign_and_CATE <- data.frame()

  for(k in 1:k_folds){

    df_learn <- df[-which(folds == k),] # train on all but kth fold

    output_learn <- learn_nuisance_k(df_learn,
                                      Y_name,
                                      A_name,
                                      sl.library.outcome,
                                      sl.library.treatment,
                                      sl.library.missingness,
                                      outcome_type,
                                      ps_trunc_level)

    k_fold_nuisance[[k]] <- output_learn[[1]]


    # if "pid" not in data, use index as pid
    if(!("pid" %in% names(df_learn))){
      df_learn$pid <- rownames(df_learn)
    }

    # note here k is actually denoting the fold that was left out of the training set
    k_fold_assign_and_CATE <- rbind(k_fold_assign_and_CATE, data.frame(pid = as.numeric(df_learn$pid),
                                                                       k = k,
                                                                       CATE_hat = output_learn[[2]]))
  }

  return(list(k_fold_nuisance, k_fold_assign_and_CATE))
}

#' Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats for kth fold
#'
#' Note training dataset contains all observations EXCEPT those assigned to kth fold
#'
#' @param df dataframe containing training dataset
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
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
learn_nuisance_k <- function(df, Y_name, A_name,
                         sl.library.outcome,
                         sl.library.treatment,
                         sl.library.missingness,
                         outcome_type,
                         ps_trunc_level = 0.01){

  # Full set of observations
  I_Y <- ifelse(is.na(df[[Y_name]]), 1, 0) #indicator for Y missing
  Y <- df[[Y_name]]
  A_vec <- df[[A_name]]
  A <- df[,A_name, drop=FALSE]
  W <- df[,!(names(df) %in% c(Y_name, A_name, "pid", "lazd90")), drop = FALSE] #QUESTION leaving laz_d90

  # Complete observations only (no NA)
  df_complete <- df[!is.na(df[[Y_name]]), ] #dataframe removing missing Ys
  Y_complete <- df_complete[[Y_name]]
  A_complete <- df_complete[, A_name, drop=FALSE]
  W_complete <- df_complete[,!(names(df) %in% c(Y_name, A_name, "pid", "lazd90")), drop = FALSE]

  ### 1 - Fit outcome model ###

  if(outcome_type == "gaussian"){
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::gaussian(),
                          cvControl = list(V=10), SL.library = sl.library.outcome)
  } else {
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::binomial(),
                          cvControl = list(V=10), SL.library = sl.library.outcome)
  }

  ### 2 - Fit treatment model ###

  pihat <- SuperLearner::SuperLearner(Y = A_vec, X = W, family = stats::binomial(),
                        cvControl = list(V=10), SL.library = sl.library.treatment)

  ### 3 - Fit missingness model ###

  deltahat <- SuperLearner::SuperLearner(Y = I_Y, X = data.frame(A, W), family = stats::binomial(),
                              cvControl = list(V=10), SL.library = sl.library.missingness)

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

  # f. Estimate CATE
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

  return(list(learned_models, CATE_hat))
}
