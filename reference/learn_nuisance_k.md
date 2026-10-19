# Estimate nuisance models (outcome, treatment, and missingness) and calculate CATE hats for kth fold

Note training dataset contains all observations EXCEPT those assigned to
kth fold

## Usage

``` r
learn_nuisance_k(
  df,
  Y_name,
  A_name,
  W_list,
  sl.library.outcome,
  sl.library.treatment,
  sl.library.missingness,
  outcome_type,
  ps_trunc_level = 0.01
)
```

## Arguments

- df:

  dataframe containing training dataset

- Y_name:

  name of outcome variable in df

- A_name:

  name of treatment variable in df

- W_list:

  character vector containing names of covariates in the dataframe to be
  used for fitting nuisance models

- sl.library.outcome:

  character vector of SuperLearner libraries to use to fit the outcome
  model

- sl.library.treatment:

  character vector of SuperLearner libraries to use to fit the treatment
  model

- sl.library.missingness:

  character vector of SuperLearner libraries to use to fit the
  missingness model

- outcome_type:

  specifying continuous (gaussian) or binomial outcome Y

- ps_trunc_level:

  numeric level to use for truncation of any predicted values that fall
  below it

## Value

- `k_fold_nuisance`:

  object of class \`Nuisance\` containing outcome, treatment, and
  missingness SuperLearner models

- `pseudo_outcome`:

  numeric vector of pseudo_outcome CATE estimates for data in df

- `shuffle_idx`:

  numeric vector corresponding to index in sorted dataframe

- `validRows`:

  SuperLearner row assignments for kth training set (training set sorted
  to have all NA outcomes first)
