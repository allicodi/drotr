# Function to learn model for CATE using k-fold cross validation

Function to learn model for CATE using k-fold cross validation

## Usage

``` r
learn_CATE(df, Z_list, k_fold_assign_and_CATE, sl.library.CATE, validRows)
```

## Arguments

- df:

  dataframe containing full dataset

- Z_list:

  character vector containing names of variables in df to use to fit
  CATE model (variables used in treatment rule)

- k_fold_assign_and_CATE:

  dataframe containing ids, fold assignments, and CATE estimate for each
  observation in df

- sl.library.CATE:

  character vector of SuperLearner libraries to use to fit the CATE
  model

- validRows:

  validRows SuperLearner row assignments from nuisance models using all
  rows in data (missingness, treatment models)

## Value

list containing SuperLearner model for CATE in each fold
