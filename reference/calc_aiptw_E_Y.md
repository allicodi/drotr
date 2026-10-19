# Function to calculate AIPTW for E_Y

Function to calculate AIPTW for E_Y

## Usage

``` r
calc_aiptw_E_Y(A, A_name, W, Y, I_Y, outcome_model, missingness_model)
```

## Arguments

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

- outcome_model:

  outcome Nuisance model

- missingness_model:

  missingness Nuisance model

## Value
