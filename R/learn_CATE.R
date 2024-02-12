#' Function to learn model for CATE using k-fold cross validation
#'
#' @param df dataframe containing full dataset
#' @param Z_list character vector containing names of variables in df to use to fit CATE model (variables used in treatment rule)
#' @param k_fold_assign_and_CATE dataframe containing ids, fold assignments, and CATE estimate for each observation in df
#' @param sl.library.CATE character vector of SuperLearner libraries to use to fit the CATE model
#' @param validRowsAll validRows SuperLearner row assignments from nuisance models using all rows in data (missingness, treatment models)
#'
#' @import SuperLearner
#' @import stats
#' @export
#'
#' @returns list containing SuperLearner model for CATE in each fold
learn_CATE <- function(df, Z_list, k_fold_assign_and_CATE, sl.library.CATE, validRowsAll){

  k_folds <- max(k_fold_assign_and_CATE$k)
  k_fold_CATE_models <- vector(mode = "list", length = k_folds)

  for(k in 1:k_folds){
    # estimate CATE hat model on all folds except kth (training data)

    # get patient ids that are NOT in the kth fold
    kth_subset_ids <- k_fold_assign_and_CATE$id[k_fold_assign_and_CATE$k == k]

    # get dataframe of patients NOT in kth fold (training data)
    df_learn <- df[(df$id %in% kth_subset_ids), , drop = FALSE]

    # get CATE hats from training set
    k_fold_assign_and_CATE_sub <- k_fold_assign_and_CATE[k_fold_assign_and_CATE$k == k,]

    # join CATEhat with data by id
    df_learn <- merge(df_learn, k_fold_assign_and_CATE_sub, by = "id")

    # df_learn$id <- as.character(df_learn$id)
    # k_fold_assign_and_CATE_sub$id <- as.character(k_fold_assign_and_CATE_sub$id)
    # df_learn <- dplyr::left_join(df_learn, k_fold_assign_and_CATE_sub, by = "id") # don't use dplyr

    validRows <- validRowsAll[[k]]

    k_fold_CATE_models[[k]] <- learn_CATE_k(df = df_learn,
                                            Z_list = Z_list,
                                            sl.library.CATE = sl.library.CATE,
                                            validRows = validRows)
  }

  return(k_fold_CATE_models)

}


#' Function to learn model for CATE in kth fold
#'
#' @param df dataframe containing training dataset
#' @param Z_list character vector containing names of variables in df to use to fit CATE model (variables used in treatment rule)
#' @param sl.library.CATE character vector of SuperLearner libraries to use to fit the CATE model
#' @param validRows SuperLearner validRows row assigments from kth fold
#'
#' @returns SuperLearner model for CATE
#'
#' @keywords internal
learn_CATE_k <- function(df, Z_list, sl.library.CATE, validRows){

  # get CATE_hat + covariates interested in using for decision rule
  Z <- df[,Z_list, drop = FALSE]
  CATE_hat <- df$pseudo_outcome

  CATE_hat_model <- SuperLearner::SuperLearner(Y = CATE_hat, X = Z, family = stats::gaussian(),
                                 cvControl = list(V=10, validRows = validRows), SL.library = sl.library.CATE)

  return(CATE_hat_model)
}
