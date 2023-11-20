#' Function to learn model for CATE using k-fold cross validation
#'
#' @param df dataframe containing full dataset
#' @param Z_list character vector containing names of variables in df to use to fit CATE model (variables used in treatment rule)
#' @param k_fold_assign_and_CATE dataframe containing pids, fold assignments, and CATE estimate for each observation in df
#' @param sl.library.CATE character vector of SuperLearner libraries to use to fit the CATE model
#'
#' @importFrom dplyr left_join
#' @import SuperLearner
#' @import stats
#' @export
#'
#' @returns SuperLearner model (discrete) for CATE
learn_CATE <- function(df, Z_list, k_fold_assign_and_CATE, sl.library.CATE){

  # if pid not in df, make pid index
  if (!"pid" %in% colnames(df)) {
    df$pid <- as.numeric(rownames(df))
  }

  k_folds <- max(k_fold_assign_and_CATE$k)
  k_fold_CATE_models <- vector(mode = "list", length = k_folds)

  for(k in 1:k_folds){
    # get patient ids that are in the kth fold
    kth_subset_pids <- k_fold_assign_and_CATE$pid[k_fold_assign_and_CATE$k == k]

    # estimate CATE hat model on all folds except kth (training data)
    df_learn <- df[!(df$pid %in% kth_subset_pids), , drop = FALSE]

    # join CATEhat with data by id
    df_learn <- dplyr::left_join(df_learn, k_fold_assign_and_CATE, by = "pid")

    k_fold_CATE_models[[k]] <- learn_CATE_k(df_learn, Z_list, sl.library.CATE)
  }

  return(k_fold_CATE_models)

}


#' Function to learn model for CATE in kth fold
#'
#' @param df dataframe containing training dataset
#' @param Z_list character vector containing names of variables in df to use to fit CATE model (variables used in treatment rule)
#' @param sl.library.CATE character vector of SuperLearner libraries to use to fit the CATE model
#'
#' @returns SuperLearner model (discrete) for CATE
#'
#' @keywords internal
learn_CATE_k <- function(df, Z_list, sl.library.CATE){

  # get CATE_hat + covariates interested in using for decision rule
  Z <- df[,Z_list, drop = FALSE]
  CATE_hat <- df$CATE_hat

  # discrete super learner final model for interpretability
  CATE_hat_model <- SuperLearner::SuperLearner(Y = CATE_hat, X = Z, family = stats::gaussian(),
                                 cvControl = list(V=10), SL.library = sl.library.CATE)

  discrete_super_learner_idx <- which.min(CATE_hat_model$cvRisk)
  CATE_hat_discrete_super_learner_model <- CATE_hat_model$fitLibrary[[discrete_super_learner_idx]]

  return(CATE_hat_discrete_super_learner_model)
}
