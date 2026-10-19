# Compare outcomes under different treatment rules

Compare outcomes under different treatment rules

## Usage

``` r
compare.otr_results(
  res_rule1,
  res_rule2,
  threshold,
  rule1_comp,
  rule2_comp,
  ...
)
```

## Arguments

- res_rule1:

  An `"otr_results"` object or `"full_otr_results"` object for treatment
  rule 1

- res_rule2:

  An `"otr_results"` object or `"full_otr_results"` object for treatment
  rule 2

- threshold:

  threshold to use for comparison of both rules (must appear in both
  otr_results objects)

- rule1_comp:

  Effect type for rule 1 to use in comparison ("treatment effect"/"te"
  or "subgroup effect"/"se")

- rule2_comp:

  Effect type for rule 2 to use in comparison ("treatment effect"/"te"
  or "subgroup effect"/"se")

- ...:

  Other arguments (not used)

## Value

dataframe containing expected value and variance for comparison
