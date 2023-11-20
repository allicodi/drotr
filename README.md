
# `drotr`

`drotr` is a package that uses doubly-robust methods to estimate optimal treatment rules (OTRs).

We define OTRs through the use of conditional average treatment effects (CATEs). The CATE can be interpreted as the expected difference in outcome under treatment vs placebo. A doubly-robust learner is used to estimate the CATE. All individuals with CATE estimates that exceed a specified threshold `t` are treated under the OTR. In other words, we only treat individuals if the treatment reduces their probability of the adverse outcome by more than `t`. 

Once observations are assigned treatment under the OTR, we can measure the risk of outcome in the optimally treated and the average treatment effects among the optimally treated. We estimate these using the Augmented Inverse Probability of Treatment Weight (AIPTW) estimators.

## Scripts:
**1. estimate_OTR.R** - This is the primary R script for finding optimal treatment rules and estimating outcomes among the optimally treated. 

**2. learn_nuisance.R** - This script is used to fit outcome, treatment, and missingness models for a given dataset. It also estimates $\hat{CATE}$ for observations in a dataset. The script can be run prior to `estimate_OTR` to pre-fit nuisance models for `estimate_OTR`, or nuisance models can be fit within `estimate_OTR`

**3. learn_CATE.R** - This script is used to fit a model for the CATE given a subset of covariates `Z` in the dataset

**4. compute_estimates.R ** -  This script is used to generate AIPTW estimators for expected outcome among the optimally treated and expected treatment effect under the optimally treated. 

## Usage:

### `estimate_EYd.R`

The `estimate_EYd.R` script contains functions necessary to estimate the expected outcome for participants in the ABCD study. The main function you'll be working with is `estimate_EYd`.

**Usage**:

Example with simulation data, pre-fit nuisance

```R
  # Example usage of `drotr` package
  # Using simulation data
  
  #
  # Simulation 1: Two normal covariates
  #
  simulate_data_1 <- function(n=1e6){
  W1 <- rnorm(n=n, 0,1)
  W2 <- rnorm(n=n, 1,1)
  
  A <- rbinom(n, 1, 0.5)
  
  Y <- W1 + W2 + W1*A + rnorm(n, 0, 1)
  
  # add missingness randomly in 10%
  delta <- rbinom(n = n, size = 1, prob = 0.1)
  
  # add missingness corresponding to delta
  Y <- ifelse(delta == 0, Y, NA)
  
  return(data.frame(Y=Y, A=A, W1=W1, W2=W2)) 
  
  }
  
  df <- simulate_data_1(n=5000)
  Z_list <- c("W1")
  sl.library <- c("SL.mean", "SL.glm", "SL.glm.interaction") #using same libraries for each step
  decision_threshold <- 0
  
  nuisance_output <- learn_nuisance(df = df,
                                    Y_name = "Y",
                                    A_name = "A",
                                    sl.library.outcome = sl.library,
                                    sl.library.treatment = sl.library,
                                    sl.library.missingness = sl.library,
                                    outcome_type = "gaussian",
                                    k_folds = 2,
                                    ps_trunc_level = 0.01)
  
  nuisance_models <- nuisance_output[[1]]
  k_fold_assign_and_CATE <- nuisance_output[[2]]
  
  results <- estimate_OTR(df = df,
                          Y_name = "Y",
                          A_name = "A",
                          Z_list = Z_list,
                          sl.library.CATE = sl.library,
                          nuisance_models = nuisance_models,
                          k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                          sl.library.outcome = sl.library,
                          sl.library.treatment = sl.library,
                          sl.library.missingness = sl.library,
                          threshold = 0,
                          k_folds = 2,
                          ps_trunc_level = 0.01,
                          outcome_type = "gaussian")


   ```
