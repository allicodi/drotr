# Make treatment decisions and compute estimated outcomes/treatment effects in kth fold

Make treatment decisions and compute estimated outcomes/treatment
effects in kth fold

## Usage

``` r
compute_estimate_k(
  df,
  Y_name,
  A_name,
  W_list,
  Z_list,
  CATE_model,
  nuisance,
  sign,
  threshold,
  ps_trunc_level = 0.01,
  truncate_CATE
)
```

## Arguments

- df:

  dataframe containing testing dataset

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

- CATE_model:

  discrete SuperLearner model for CATE

- nuisance:

  object of class \`Nuisance\` containing outcome, treatment, and
  missingness SuperLearner models

- sign:

  sign of threshold, + indicates \`Y\` is a desirable outcome, -
  indicates \`Y\` is an undesirable outcome

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
