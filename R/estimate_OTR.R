#' Main function to calculate estimated treatment effects for treatment rule
#'
#' Primary function to estimate treatment effects for a given treatment rule.
#' Can fit nuisance models internally or be provided with pre-fit nuisance models for the given dataset.
#'
#' @param df dataframe containing full dataset
#' @param Y_name name of outcome variable in df
#' @param A_name name of treatment variable in df
#' @param W_list character vector containing names of covariates in the dataframe to be used for fitting nuisance models
#' @param Z_list character vector containing names of variables in df used to fit CATE model (variables used in treatment rule)
#' @param id_name name of participant ID variable
#' @param sl.library.CATE character vector of SuperLearner libraries to use to fit the CATE model
#' @param nuisance_models list of objects of class `Nuisance` containing outcome, treatment, and missingness SuperLearner models (only include if using pre-fit nuisance models)
#' @param k_fold_assign_and_CATE dataframe containing ids, fold assignments, and CATE estimate for each observation in df (only include if using pre-fit nuisance models)
#' @param validRows list containing validRows assignents for nuisance models using all observations (only include if using pre-fit nuisance models)
#' @param sl.library.outcome character vector of SuperLearner libraries to use to fit the outcome models
#' @param sl.library.treatment character vector of SuperLearner libraries to use to fit the treatment models
#' @param sl.library.missingness character vector of SuperLearner libraries to use to fit the missingness models
#' @param threshold character vector of decision thresholds for CATE to determine OTR. Values should be positive if `Y_name` is desirable outcome, negative if `Y_name` is undesirable outcome. If threshold is 0, use +0 for desirable, -0 for undesirable.
#' @param k_folds integer number of folds to use for cross-validation (must specify if fitting outcome, treatment, and missingness models. Otherwise uses k from `k_fold_assign_and_CATE`)
#' @param ps_trunc_level numeric evel below which propensity scores will be truncated (to avoid errors in computing AIPTW)
#' @param outcome_type outcome_type specifying continuous (outcome_type = "gaussian") or binary (outcome_type = "binomial") outcome Y (if not providing pre-fit nuisance models)
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
estimate_OTR <- function(df,
                         Y_name,
                         A_name,
                         W_list,
                         Z_list,
                         id_name = NULL,
                         sl.library.CATE,
                         nuisance_models = NULL,
                         k_fold_assign_and_CATE = NULL,
                         validRows = NULL,
                         sl.library.outcome = NULL,
                         sl.library.treatment = NULL,
                         sl.library.missingness = NULL,
                         threshold = c("0.05"),
                         k_folds = 2,
                         ps_trunc_level = 0.01,
                         outcome_type = "gaussian"){

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
  # 1 - Fit nuisance models (if not provided)
  # --------------------------------------------------------------------------

  # Check that both models and fold assignments are provided, if one provided, exit
  if ((is.null(nuisance_models) & !is.null(k_fold_assign_and_CATE)) | (!is.null(nuisance_models) & is.null(k_fold_assign_and_CATE))){
    return(print("Must provide both nuisance models and fold assignments to use for estimation."))
  }
  else if (is.null(nuisance_models) & is.null(k_fold_assign_and_CATE)) { # If neither provided, estimate models & assign folds

    # Verify all necessary SuperLearner libraries were input
    if(is.null(sl.library.outcome) | is.null(sl.library.treatment) | is.null(sl.library.missingness)){
      return(print("Must provide outcome, treatment, and missingness libraries to estimate nuisance models."))
    }

    nuisance_output <- drotr::learn_nuisance(df, Y_name, A_name, W_list, id_name, sl.library.outcome, sl.library.treatment,
                                      sl.library.missingness, outcome_type, k_folds, ps_trunc_level)

    nuisance_models <- nuisance_output$nuisance_models
    k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
    validRows <- nuisance_output$validRows

  }

  # --------------------------------------------------------------------------
  # 2 - Learn CATE models
  # --------------------------------------------------------------------------
  CATE_models <- drotr::learn_CATE(df, Z_list, k_fold_assign_and_CATE, sl.library.CATE, validRows)

  # --------------------------------------------------------------------------
  # 3 - Make treatment decisions and compute treatment effects
  # --------------------------------------------------------------------------
  results <- drotr::compute_estimates(df, Y_name, A_name, W_list, Z_list,
                                      k_fold_assign_and_CATE,
                                      nuisance_models, CATE_models,
                                      threshold, ps_trunc_level)

  results <- list(results)
  names(results) <- "results"

  results$nuisance_models <- nuisance_models
  results$CATE_models <- CATE_models
  results$Z_list <- Z_list

  class(results) <- "otr_results"

  return(results)

}
