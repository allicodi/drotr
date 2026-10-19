# Helper function to average results across multiple seeds

Helper function to average results across multiple seeds

## Usage

``` r
average_across_seeds(results_list, threshold)
```

## Arguments

- results_list:

  list of results using same models & Z_list, different seeds

- threshold:

  threshold of results to use for comparison

## Value

summary_df dataframe summarizing results across length(results_list)
seeds
