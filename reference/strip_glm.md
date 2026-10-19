# Helper function to remove unnecessary output from SuperLearner model to reduce output size

Reference:
https://www.r-bloggers.com/2019/12/need-to-save-rs-lm-or-glm-models-trim-the-fat/

## Usage

``` r
strip_glm(cm, earth = FALSE)
```

## Arguments

- cm:

  SuperLearner GLM or Earth model

- earth:

  flag for if model is SuperLearner Earth model, defaults to GLM

## Value

model with unnecessary features removed
