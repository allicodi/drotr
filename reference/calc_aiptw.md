# Function to calculate AIPTW

Function to calculate AIPTW

## Usage

``` r
calc_aiptw(
  a,
  A,
  A_name,
  W,
  Y,
  I_Y,
  d_pred,
  mean_dZ,
  outcome_model,
  treatment_model,
  missingness_model,
  ps_trunc_level,
  idx_sub
)
```

## Arguments

- a:

  value of treatment to set for all observations

- A:

  vector of observed treatments

- A_name:

  name of treatment variable

- W:

  dataframe of covariates

- Y:

  vector of observed outcome variable

- I_Y:

  indicator that Y is missing (I_Y = 1 if missing, I_Y = 0 if observed)

- d_pred:

  vector of treatment decisions (d_pred = 1 if recommend treatment under
  rule)

- mean_dZ:

  proportion recommended treatment under rule

- outcome_model:

  outcome Nuisance model

- treatment_model:

  treatment Nuisance model

- missingness_model:

  missingness Nuisance model

- ps_trunc_level:

  threshold to truncate propensity score

- idx_sub:

  indices of observations that are recommended treatment (d_pred == 1)

## Value
