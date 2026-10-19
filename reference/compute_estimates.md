# Make treatment decisions and compute estimated outcomes/treatment effects

Make treatment decisions and compute estimated outcomes/treatment
effects

## Usage

``` r
compute_estimates(
  df,
  Y_name,
  A_name,
  W_list,
  Z_list,
  k_fold_assign_and_CATE,
  nuisance_models,
  CATE_models,
  threshold,
  ps_trunc_level = 0.01,
  truncate_CATE
)
```

## Arguments

- df:

  dataframe containing dataset

- Y_name:

  name of outcome variable in df

- A_name:

  name of treatment variable in df

- W_list:

  character vector containing names of covariates in the dataframe to be
  used in nuisance models

- Z_list:

  character vector containing names of variables in df used to fit CATE
  model (variables used in treatment rule)

- k_fold_assign_and_CATE:

  dataframe containing ids, fold assignments, and CATE estimate for each
  observation in df

- nuisance_models:

  list of objects of class \`Nuisance\` containing outcome, treatment,
  and missingness SuperLearner models

- CATE_models:

  list of discrete SuperLearner models for CATE from each fold

- threshold:

  character vector of decision thresholds for CATE to determine OTR.
  Values should be positive if \`Y_name\` is desirable outcome, negative
  if \`Y_name\` is undesirable outcome. If threshold is 0, use +0 for
  desirable, -0 for undesirable.

- ps_trunc_level:

  numeric level to use for truncation of any predicted values that fall
  below it

- truncate_CATE:

  logical to indicate if large CATE predictions should be truncated at
  -1 and 1 (default = TRUE)

## Value

- `aggregated_results`:

  dataframe of overall results aggregated across \`k\` folds for given
  threshold

- `k_fold_results`:

  list of results by fold for given threshold

- `decision_df`:

  original dataset with decision made for each observation at a given
  threshold

- `k_non_na`:

  folds that did not have proportion treated = 1 or = 0 (causing some
  effect estimates to be NA)
