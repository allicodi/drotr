# Function to apply existing OTR(s) to external data

Function to apply existing OTR(s) to external data

## Usage

``` r
apply_OTR(
  df,
  CATE_models = NULL,
  d_pred = NULL,
  Y_name,
  A_name,
  W_list,
  Z_list = NULL,
  id_name = NULL,
  threshold = c("0.05"),
  nuisance_models = NULL,
  sl.library.outcome = NULL,
  sl.library.treatment = NULL,
  sl.library.missingness = NULL,
  k_folds = 5,
  ps_trunc_level = 0.01,
  outcome_type = "gaussian",
  truncate_CATE = "TRUE"
)
```

## Arguments

- df:

  dataframe containing external dataset to apply rule(s) to; must
  contain Z_list variables that are same as pre-trained rule(s)

- CATE_models:

  list of CATE model(s) to apply to external dataset (if applying
  existing rule to new data)

- d_pred:

  vector of binary treatment decisions (if have existing decisions and
  just want to estimate effects)

- Y_name:

  name of outcome variable. Outcome variable in rule being applied
  should be the same as outcome variable in new data.

- A_name:

  name of treatment variable. Treatment variable in rule being applied
  should be the same as the treatment variable in the new data.

- W_list:

  character vector containing names of covariates in the dataframe used
  for nuisance models.

- Z_list:

  character vector containing names of variables in df used to fit CATE
  model (variables used in treatment rule; must be same names as used in
  pre-fit CATE model(s); only required if applying rule)

- id_name:

  name of participant ID variable if present in data

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
