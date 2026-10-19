# Function to compute estimates on external dataset

Function to compute estimates on external dataset

## Usage

``` r
compute_estimates_external(
  df,
  Y_name,
  A_name,
  W_list,
  Z_list,
  CATE_preds = NULL,
  d_pred = NULL,
  threshold = 0,
  nuisance_models = NULL,
  sl.library.outcome = NULL,
  sl.library.treatment = NULL,
  sl.library.missingness = NULL,
  k_folds = 5,
  ps_trunc_level = 0.01,
  outcome_type = "gaussian",
  id_name = NULL
)
```

## Arguments

- df:

  dataframe containing external dataset to apply rule(s) to; must
  contain Z_list variables that are same as pre-trained rule(s)

- Y_name:

  name of outcome variable in df

- A_name:

  name of treatment variable in df

- W_list:

  character vector containing names of covariates in the dataframe to be
  used for fitting nuisance models

- Z_list:

  character vector containing names of variables in df used to fit CATE
  model (variables used in treatment rule; must be same names as used in
  pre-fit CATE model(s))

- CATE_preds:

  dataframe with CATE predictions for each CATE model & truncation flags
  for extreme predictions

- d_pred:

  vector of binary treatment decisions (if have existing decisions and
  just want to estimate effects)

- threshold:

  character vector of decision thresholds for CATE to determine OTR.
  Values should be positive if \`Y_name\` is desirable outcome, negative
  if \`Y_name\` is undesirable outcome. If threshold is 0, use +0 for
  desirable, -0 for undesirable.

- nuisance_models:

  list of objects of class \`Nuisance\` containing outcome, treatment,
  and missingness SuperLearner models (only include if using pre-fit
  nuisance models)

- sl.library.outcome:

  character vector of SuperLearner libraries to use to fit the outcome
  models (if fitting nuisance internally)

- sl.library.treatment:

  character vector of SuperLearner libraries to use to fit the treatment
  models (if fitting nuisance internally)

- sl.library.missingness:

  character vector of SuperLearner libraries to use to fit the
  missingness models (if fitting nuisance internally)

- k_folds:

  integer number of folds to use for nuisance model cross-validation
  (must specify if fitting nuisance models here)

- ps_trunc_level:

  numeric level below which propensity scores will be truncated (to
  avoid errors in computing AIPTW)

- outcome_type:

  outcome_type specifying continuous (outcome_type = "gaussian") or
  binary (outcome_type = "binomial") outcome Y (if not providing pre-fit
  nuisance models)

- id_name:

  name of participant ID variable if present in data

## Value

- `aggregated_results`:

  Results averaged over each nuisance model

- `k_fold_results`:

  Results for each nuisance model

- `decision_df`:

  Dataframe with CATE prediction from every individual CATE_model, flag
  for CATE truncation, average CATE predition used for treatment
  decision, treatment decision

- `k_non_na`:

  equal to number of k_folds used for nuisance models. Ignore–
  placeholder to keep same format as Results object from internal
  compute estimates
