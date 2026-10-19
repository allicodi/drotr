# Main function to calculate estimated treatment effects for treatment rule

Primary function to estimate treatment effects for a given treatment
rule. Can fit nuisance models internally or be provided with pre-fit
nuisance models for the given dataset.

## Usage

``` r
estimate_OTR(
  df,
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
  outcome_type = "gaussian",
  truncate_CATE = "TRUE"
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

- Z_list:

  character vector containing names of variables in df used to fit CATE
  model (variables used in treatment rule)

- id_name:

  name of participant ID variable

- sl.library.CATE:

  character vector of SuperLearner libraries to use to fit the CATE
  model

- nuisance_models:

  list of objects of class \`Nuisance\` containing outcome, treatment,
  and missingness SuperLearner models (only include if using pre-fit
  nuisance models)

- k_fold_assign_and_CATE:

  dataframe containing ids, fold assignments, and CATE estimate for each
  observation in df (only include if using pre-fit nuisance models)

- validRows:

  list containing validRows assignents for nuisance models using all
  observations (only include if using pre-fit nuisance models)

- sl.library.outcome:

  character vector of SuperLearner libraries to use to fit the outcome
  models

- sl.library.treatment:

  character vector of SuperLearner libraries to use to fit the treatment
  models

- sl.library.missingness:

  character vector of SuperLearner libraries to use to fit the
  missingness models

- threshold:

  character vector of decision thresholds for CATE to determine OTR.
  Values should be positive if \`Y_name\` is desirable outcome, negative
  if \`Y_name\` is undesirable outcome. If threshold is 0, use +0 for
  desirable, -0 for undesirable.

- k_folds:

  integer number of folds to use for cross-validation (must specify if
  fitting outcome, treatment, and missingness models. Otherwise uses k
  from \`k_fold_assign_and_CATE\`)

- ps_trunc_level:

  numeric evel below which propensity scores will be truncated (to avoid
  errors in computing AIPTW)

- outcome_type:

  outcome_type specifying continuous (outcome_type = "gaussian") or
  binary (outcome_type = "binomial") outcome Y (if not providing pre-fit
  nuisance models)

- truncate_CATE:

  logical to indicate if large CATE predictions should be truncated at
  -1 and 1 (default = TRUE)

## Value

- `results`:

  list of \`Results\` objects for each threshold. See description of
  \`Results\` object in \`compute_estimates\`

- `nuisance_models`:

  list of \`Nuisance\` objects containing outcome, treatment, and
  missingness models used in each fold

- `CATE_models`:

  CATE model used in each fold

- `Z_list`:

  character vector containing names of variables in df used to fit CATE
  model (variables used in treatment rule)
