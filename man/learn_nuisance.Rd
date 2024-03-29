% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/learn_nuisance.R
\name{learn_nuisance}
\alias{learn_nuisance}
\title{Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats using k-fold cross validation}
\usage{
learn_nuisance(
  df,
  Y_name,
  A_name,
  W_list,
  id_name = NULL,
  sl.library.outcome,
  sl.library.treatment,
  sl.library.missingness,
  outcome_type,
  k_folds = 2,
  ps_trunc_level = 0.01
)
}
\arguments{
\item{df}{dataframe containing full dataset}

\item{Y_name}{name of outcome variable in df}

\item{A_name}{name of treatment variable in df}

\item{W_list}{character vector containing names of covariates in the dataframe to be used for fitting nuisance models}

\item{id_name}{name of patient id variable in dataset if applicable, will default to NULL and use observation index}

\item{sl.library.outcome}{character vector of SuperLearner libraries to use to fit the outcome models}

\item{sl.library.treatment}{character vector of SuperLearner libraries to use to fit the treatment models}

\item{sl.library.missingness}{character vector of SuperLearner libraries to use to fit the missingness models}

\item{outcome_type}{specifying continuous (outcome_type = "gaussian") or binary (outcome_type = "binomial") outcome Y}

\item{k_folds}{number of folds for k_fold cross validation}

\item{ps_trunc_level}{numeric level to use for truncation of any predicted values that fall below it}
}
\value{
\describe{
 \item{\code{k_fold_nuisance}}{list of `Nuisance` objects (fit nuisance models) for each of k folds}
 \item{\code{k_fold_assign_and_CATE}}{dataframe of CATE estimates, k-1 folds, pseudo-outcome, and shuffle idx corresponding to validRows for each observation}
 \item{\code{validRows}}{list of innerCV SuperLearner row assignments for each training set}
 \item{\code{fold_assignments}}{dataframe containing fold assignments for each id}
 }
}
\description{
Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats using k-fold cross validation
}
