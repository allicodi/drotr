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
#'  \item{\code{k_fold_assign_and_CATE}}{dataframe of CATE estimates, k-1 folds, pseudo-outcome, and shuffle idx corresponding to validRows for each observation}
#'  \item{\code{validRows}}{list of innerCV SuperLearner row assignments for each training set}
#'  \item{\code{fold_assignments}}{dataframe containing fold assignments for each id}
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
    df$id <- 1:nrow(df)
  } else {
    df$id <- df[[id_name]]
  }

  # split into k equal sized folds and randomly shuffle
  folds <- cut(seq(1,nrow(df)),breaks=k_folds,labels=FALSE)
  folds <- sample(folds)

  fold_assignments <- data.frame(id = df$id, fold = folds)

  k_fold_nuisance <- vector(mode = "list", length = k_folds)
  k_fold_assign_and_CATE <- data.frame()

  validRowsByFold <- vector(mode = "list", length = k_folds)

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
                                                                       pseudo_outcome = output_learn$pseudo_outcome,
                                                                       shuffle_idx = output_learn$shuffle_idx))
    validRowsByFold[[k]] <- output_learn$validRows
  }

  return(list(nuisance_models = k_fold_nuisance,
              k_fold_assign_and_CATE = k_fold_assign_and_CATE,
              validRows = validRowsByFold,
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
#' @param outcome_type specifying continuous (gaussian) or binomial outcome Y
#' @param ps_trunc_level numeric level to use for truncation of any predicted values that fall below it
#'
#' @return
#' \describe{
#'  \item{\code{k_fold_nuisance}}{object of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models}
#'  \item{\code{pseudo_outcome}}{numeric vector of pseudo_outcome CATE estimates for data in df}
#'  \item{\code{shuffle_idx}}{numeric vector corresponding to index in sorted dataframe}
#'  \item{\code{validRows}}{SuperLearner row assignments for kth training set (training set sorted to have all NA outcomes first)}
#'  }
#' @keywords internal
learn_nuisance_k <- function(df, Y_name, A_name, W_list,
                         sl.library.outcome,
                         sl.library.treatment,
                         sl.library.missingness,
                         outcome_type,
                         ps_trunc_level = 0.01){

  # assign var original indices
  df$original_id <- 1:nrow(df)

  # get indices in original dataframe corresponding to NA, non-NA
  idx_na_Y <- which(is.na(df[[Y_name]]))
  idx_no_na_Y <- which(!is.na(df[[Y_name]]))

  # get number of NA observations, overall observations
  nNA <- length(idx_na_Y)
  n <- nrow(df)

  # sort dataframe to put all NAs at beginning
  df_sort <- df[c(idx_na_Y, idx_no_na_Y), ]

  # assign var sorted indices
  df_sort$df_sort_id <- 1:nrow(df_sort)

  # Full set of observations
  I_Y <- ifelse(is.na(df_sort[[Y_name]]), 1, 0) #indicator for Y missing
  Y <- df_sort[[Y_name]]
  A_vec <- df_sort[[A_name]]
  A <- df_sort[,A_name, drop=FALSE]
  W <- df_sort[,W_list, drop=FALSE]

  # Complete observations only (no NA)
  df_complete <- df_sort[!is.na(df_sort[[Y_name]]), ] #dataframe removing missing Ys
  Y_complete <- df_complete[[Y_name]]
  A_complete <- df_complete[, A_name, drop=FALSE]
  W_complete <- df_complete[, W_list, drop = FALSE]

  # split NA observations into validRows folds (if NA present)
  # changing V from default 10 --> 3 to have more observations for fitting CATE model later on
  if(nNA > 0){
    rows_na_Y <- split(sample(1:nNA), rep(1:3, length = nNA))
  }else{
    rows_na_Y <- vector(mode = "list", length = 0L)
  }

  # split remaining observations into validRows folds
  rows_no_na_Y <- split(sample((nNA+1):n), rep(1:3, length = n-nNA))
  master_validRows <- vector(mode = "list", length = 3)
  for (vv in seq_len(3)) {
    if(length(rows_na_Y) >= vv & length(rows_no_na_Y) >= vv){
      master_validRows[[vv]] <- c(rows_na_Y[[vv]], rows_no_na_Y[[vv]])
    }else if(length(rows_na_Y) < vv){
      master_validRows[[vv]] <- rows_no_na_Y[[vv]]
    }else if(length(rows_no_na_Y) < vv){
      master_validRows[[vv]] <- rows_na_Y[[vv]]
    }
  }

  ### 1 - Fit outcome model ###

  # muhat fit using only non-NA observations (use validRows excluding NAs)
  muhat_validRows <- lapply(master_validRows, function(x){
    col <- x - nNA
    col[col > 0]
  })

  if(outcome_type == "gaussian"){
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::gaussian(),
                          cvControl = list(validRows = muhat_validRows, V = 3), SL.library = sl.library.outcome, control = list(saveCVFitLibrary = TRUE))
  } else {
    muhat <- SuperLearner::SuperLearner(Y = Y_complete, X = data.frame(A_complete, W_complete), family = stats::binomial(),
                          cvControl = list(validRows = muhat_validRows, V = 3), SL.library = sl.library.outcome, control = list(saveCVFitLibrary = TRUE))
  }

  muhat.cvFitLibrary <- muhat$cvFitLibrary
  muhat.coef <- muhat$coef

  ### 2 - Fit treatment model ###

  pihat <- SuperLearner::SuperLearner(Y = A_vec, X = W, family = stats::binomial(),
                        cvControl = list(validRows = master_validRows, V = 3), SL.library = sl.library.treatment, control = list(saveCVFitLibrary = TRUE))

  pihat.cvFitLibrary <- pihat$cvFitLibrary
  pihat.coef <- pihat$coef

  ### 3 - Fit missingness model ###

  deltahat <- SuperLearner::SuperLearner(Y = I_Y, X = data.frame(A, W), family = stats::binomial(),
                              cvControl = list(validRows = master_validRows, V = 3), SL.library = sl.library.missingness, control = list(saveCVFitLibrary = TRUE))

  deltahat.cvFitLibrary <- deltahat$cvFitLibrary
  deltahat.coef <- deltahat$coef

  ### 4 - Create pseudo-outcome  ###

  n_muhat_lib <- length(muhat.cvFitLibrary[[1]])
  n_pihat_lib <- length(pihat.cvFitLibrary[[1]])
  n_deltahat_lib <- length(deltahat.cvFitLibrary[[1]])

  muhat_obs_matrix <- matrix(NA, nrow = nrow(df), ncol = n_muhat_lib)
  muhat_1_matrix <- matrix(NA, nrow = nrow(df), ncol = n_muhat_lib)
  muhat_0_matrix <- matrix(NA, nrow = nrow(df), ncol = n_muhat_lib)
  pihat_matrix <- matrix(NA, nrow = nrow(df), ncol = n_pihat_lib)
  deltahat_matrix <- matrix(NA, nrow = nrow(df), ncol = n_deltahat_lib)

  for(v in 1:3){
    muhat.v <- muhat.cvFitLibrary[[v]]
    pihat.v <- pihat.cvFitLibrary[[v]]
    deltahat.v <- deltahat.cvFitLibrary[[v]]

    # a. get pred from outcome model under obs trt and obs cov (outcome model)
    for(model in 1:length(muhat.v)){
      muhat.pred <- stats::predict(
        muhat.v[[model]],
        newdata = data.frame(A, W)[master_validRows[[v]], , drop = FALSE],
        family = muhat$family
      )
      muhat_obs_matrix[master_validRows[[v]], model] <- muhat.pred
    }

    # b. predict from trt = 1 and obs cov (outcome model)
    om_trt1_df <- data.frame(1, W)
    names(om_trt1_df)[1] <- A_name

    for(model in 1:length(muhat.v)){
      muhat1.pred <- stats::predict(
        muhat.v[[model]],
        newdata = om_trt1_df[master_validRows[[v]], , drop = FALSE],
        family = muhat$family
      )
      muhat_1_matrix[master_validRows[[v]], model] <- muhat1.pred
    }

    # c. predict from trt = 0 and obs cov (outcome model)
    om_trt0_df <- data.frame(0, W)
    names(om_trt0_df)[1] <- A_name

    for(model in 1:length(muhat.v)){
      muhat0.pred <- stats::predict(
        muhat.v[[model]],
        newdata = om_trt0_df[master_validRows[[v]], , drop = FALSE],
        family = muhat$family
      )
      muhat_0_matrix[master_validRows[[v]], model] <- muhat0.pred
    }

    # d. predict treatment from obs cov & truncate observations that are too small (treatment model)
    for(model in 1:length(pihat.v)){
      pihat.pred <- stats::predict(
        pihat.v[[model]],
        newdata = data.frame(W)[master_validRows[[v]], , drop = FALSE],
        family = pihat$family
      )

      # truncate if probability too close to 0/1
      pihat.pred[pihat.pred < ps_trunc_level] <- ps_trunc_level
      pihat.pred[pihat.pred > 1 - ps_trunc_level] <- 1 - ps_trunc_level

      pihat_matrix[master_validRows[[v]], model] <- pihat.pred
    }

    # e. predict missingness under obs cov, obs trt (missingness model)
    for(model in 1:length(deltahat.v)){
      deltahat.pred <- stats::predict(
        deltahat.v[[model]],
        newdata = data.frame(A, W)[master_validRows[[v]], , drop = FALSE],
        family = deltahat$family
      )

      # truncate if probability too close to 0/1
      deltahat.pred[deltahat.pred < ps_trunc_level] <- ps_trunc_level

      deltahat_matrix[master_validRows[[v]], model] <- deltahat.pred
    }

  }

  muhat_obs_weighted <- muhat_obs_matrix %*% muhat.coef
  muhat_1_weighted <- muhat_1_matrix %*% muhat.coef
  muhat_0_weighted <- muhat_0_matrix %*% muhat.coef
  pihat_weighted <- pihat_matrix %*% pihat.coef
  deltahat_weighted <- deltahat_matrix %*% deltahat.coef

  # check that this is the same as passing in correct models
  sl.outcome.correct <- function(X, an_grp_01){
    -0.4213183637 +
      0.0079037187 * X$dy1_scrn_diardays +
      0.1150761291 * I(X$site == "Kenya") +
      -0.0539769251 * I(X$site == "Malawi") +
      0.0788648620 * I(X$site == "Mali") +
      -0.0250071035 * I(X$site == "India") +
      -0.1211349999 * I(X$site == "Tanzania") +
      -0.0058588270 * I(X$site == "Pakistan") +
      -0.0626945585 * I(X$an_ses_quintile == "2nd quintile of SES") +
      -0.0246203478 * I(X$an_ses_quintile == "3rd quintile of SES") +
      0.0214125476 * I(X$an_ses_quintile == "4th quintile of SES") +
      0.0517797108  * I(X$an_ses_quintile == "5th quintile of SES") +
      -0.0001681260 * X$agemchild +
      0.8384232165 * X$lfazscore +
      -0.0763975679 * I(X$shigella_new > 0) +
      0.1940622 * I(I(X$shigella_new > 0) * an_grp_01) +
      -0.007401625 * I(I(X$shigella_new > 0) * X$lfazscore * an_grp_01) +
      -0.0025 * I(I(X$shigella_new > 0) * X$lfazscore * X$agemchild * an_grp_01)
  }

  # test hardcode
  # muhat_obs_weighted <- sl.outcome.correct(X = df_sort, an_grp_01 = df_sort$an_grp_01)
  # muhat_1_weighted <- sl.outcome.correct(X = df_sort, an_grp_01 = 1)
  # muhat_0_weighted <- sl.outcome.correct(X = df_sort, an_grp_01 = 0)
  # pihat_weighted <- rep(0.5, nrow(df_sort))
  # deltahat_weighted <- rep(0.04, nrow(df_sort))

  # f. Estimate CATE (pseudo-outcome)
  p1 <- ((2*A_vec - 1)*as.numeric(!I_Y) ) / (pihat_weighted * (1-deltahat_weighted))
  p2 <- Y - muhat_obs_weighted
  p2 <- ifelse(is.na(p2), 0, p2)
  p3 <- muhat_1_weighted - muhat_0_weighted

  # once confirmed here, remove hard coding, use wrappers for increasing sample size on the cluster
  CATE_hat <- p1*p2 + p3 + rnorm(nrow(df_sort), 0, 1/sqrt(nrow(df_sort)))

  # reorder CATE hat to correspond to original dataframe
  df_sort$CATE_hat <- CATE_hat
  reorder_df <- df_sort[order(df_sort$original_id),]

  # eliminate parts of glm that are too large
  muhat <- strip_nuisance(muhat)
  pihat <- strip_nuisance(pihat)
  deltahat <- strip_nuisance(deltahat)

  # Create Nuisance object
  learned_models <- list(outcome_model = muhat,
                         treatment_model = pihat,
                         missingness_model = deltahat)

  class(learned_models) <- "Nuisance"

  return(list(nuisance_models = learned_models,
              pseudo_outcome = reorder_df$CATE_hat,
              shuffle_idx = reorder_df$df_sort_id,
              validRows = master_validRows))
}
