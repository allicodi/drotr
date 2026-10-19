# Function to learn model for CATE in kth fold

Function to learn model for CATE in kth fold

## Usage

``` r
learn_CATE_k(df, Z_list, sl.library.CATE, validRows)
```

## Arguments

- df:

  dataframe containing training dataset

- Z_list:

  character vector containing names of variables in df to use to fit
  CATE model (variables used in treatment rule)

- sl.library.CATE:

  character vector of SuperLearner libraries to use to fit the CATE
  model

- validRows:

  SuperLearner validRows row assigments from kth fold

## Value

SuperLearner model for CATE
