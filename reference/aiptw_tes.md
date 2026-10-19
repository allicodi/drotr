# Function to estimate treatment effects using AIPTW for external data

Function to estimate treatment effects using AIPTW for external data

## Usage

``` r
aiptw_tes(
  df,
  d_pred,
  Y_name,
  A_name,
  W_list,
  Z_list,
  nuisance_models,
  ps_trunc_level
)
```

## Arguments

- df:

  dataframe containing external dataset to apply rule(s) to; must
  contain Z_list variables that are same as pre-trained rule(s)

- d_pred:

  treatment decisions from single CATE model and threshold combo

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

- nuisance_models:

  list of objects of class \`Nuisance\` containing outcome, treatment,
  and missingness SuperLearner models

- ps_trunc_level:

  numeric evel below which propensity scores will be truncated (to avoid
  errors in computing AIPTW)

## Value

list of effect estimates, influence function matrix
