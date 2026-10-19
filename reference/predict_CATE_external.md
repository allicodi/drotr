# Function to predict CATEs & get treatment decisions on external dataset

Function to predict CATEs & get treatment decisions on external dataset

## Usage

``` r
predict_CATE_external(df, CATE_models, Z_list, truncate_CATE = TRUE)
```

## Arguments

- df:

  dataframe containing external dataset to apply rule(s) to; must
  contain Z_list variables that are same as pre-trained rule(s)

- CATE_models:

  list of CATE model(s) to apply to external dataset

- Z_list:

  character vector containing names of variables in df used to fit CATE
  model (variables used in treatment rule; must be same names as used in
  pre-fit CATE model(s))

- truncate_CATE:

  logical to indicate if large CATE predictions should be truncated at
  -1 and 1 (default = TRUE)

## Value

dataframe with columns for id, CATEs for each model, and average CATE
for rule
