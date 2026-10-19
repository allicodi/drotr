# Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats using k-fold cross validation

Estimate nuisance models (outcome, treatment, and missingness) and
calculate CATE hats using k-fold cross validation

## Usage

``` r
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
```

## Arguments

- df:

  dataframe containing full dataset

- Y_name:

  name of outcome variable in df

- A_name:

  name of treatment variable in df

- W_list:

  character vector containing names of covariates in the dataframe to be
  used for fitting nuisance models

- id_name:

  name of patient id variable in dataset if applicable, will default to
  NULL and use observation index

- sl.library.outcome:

  character vector of SuperLearner libraries to use to fit the outcome
  models

- sl.library.treatment:

  character vector of SuperLearner libraries to use to fit the treatment
  models

- sl.library.missingness:

  character vector of SuperLearner libraries to use to fit the
  missingness models

- outcome_type:

  specifying continuous (outcome_type = "gaussian") or binary
  (outcome_type = "binomial") outcome Y

- k_folds:

  number of folds for k_fold cross validation

- ps_trunc_level:

  numeric level to use for truncation of any predicted values that fall
  below it

## Value

- `k_fold_nuisance`:

  list of \`Nuisance\` objects (fit nuisance models) for each of k folds

- `k_fold_assign_and_CATE`:

  dataframe of CATE estimates, k-1 folds, pseudo-outcome, and shuffle
  idx corresponding to validRows for each observation

- `validRows`:

  list of innerCV SuperLearner row assignments for each training set

- `fold_assignments`:

  dataframe containing fold assignments for each id
