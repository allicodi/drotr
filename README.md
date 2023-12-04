
# `drotr`

`drotr` is a package that uses doubly-robust methods to estimate optimal treatment rules (OTRs).

We define OTRs through the use of conditional average treatment effects (CATEs). The CATE can be interpreted as the expected difference in outcome under treatment vs placebo. A doubly-robust learner is used to estimate the CATE. All individuals with CATE estimates that exceed a specified threshold 't' are treated under the OTR. In other words, we only treat individuals if the treatment reduces their probability of the adverse outcome by more than 't'. 

Once observations are assigned treatment under the OTR, we can measure the risk of outcome in the optimally treated and the average treatment effects among the optimally treated. We estimate these using Augmented Inverse Probability of Treatment Weight (AIPTW) estimators.

## Installation:

A developmental release may be installed from GitHub via devtools with:

```devtools::install_github("allicodi/drotr")```

## Usage:

Suppose we have a dataset `df` consisting of continuous or binary outcome variable `Y`, baseline covariates `W`, and binary treatment variable `A`. We want to estimate an optimal treatment rule (OTR) based on a subest of covariates `Z`. `Y` may also contain missing data. 

```R
  
  #
  # Simulation with two normal covariates
  #
  simulate_data <- function(n=1e6){
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
  
  df <- simulate_data(n=5000)
  Z_list <- c("W1")
  
```

The function `estimate_OTR` will assign treatment to all observations in `df` with a clinically relevant individual level treatment effect. It will also estimate the treatment effect of those assigned treatment by the OTR. 

```R
  
  # Nuisance model SuperLearner libraries
  sl.library.outcome <- c("SL.mean", "SL.glm", "SL.glm.interaction")      # libraries to use for outcome model
  sl.library.treatment <- c("SL.mean", "SL.glm", "SL.glm.interaction")    # libraries to use for treatment model
  sl.library.missingness <- c("SL.mean", "SL.glm", "SL.glm.interaction")  # libraries to use for missingness model
  
  # CATE model SuperLearner libraries
  sl.library.CATE <- c("SL.mean", "SL.glm", "SL.glm.interaction")
  
  # List of covariates to use to estimate Nuisance models
  W_list <- c("W1", "W2")
  
  # List of covariates to use to estimate CATE model 
  Z_list <- c("W1")
  
  decision_threshold <- 0  # clinically relevance threshold for treatment effect (>=0 if desired outcome Y, negative in undesirable)
  
  results <- estimate_OTR(df = df,
                          Y_name = "Y",
                          A_name = "A",
                          W_list = W_list,
                          Z_list = Z_list,
                          sl.library.CATE = sl.library.CATE,
                          sl.library.outcome = sl.library.outcome,
                          sl.library.treatment = sl.library.treatment,
                          sl.library.missingness = sl.library.missingness,
                          threshold = 0,
                          k_folds = 2,
                          ps_trunc_level = 0.01,
                          outcome_type = "gaussian")
                          
    overall_results <- results$overall_results      # dataframe of overall results aggregated across `k` folds
    EY_A1_d1 <- results$EY_A1_d1                    # dataframe of AIPTW for optimally treated in each fold
    EY_A0_d1 <- results$EY_A0_d1                    # dataframe of AIPTW for not treating those who should be treated under decision rule in each fold
    treatment_effect <- results$treatment_effect    # dataframe of treatment effect in each fold
    decision_df <- results$decision_df              # original dataset with decision made for each observation
    CATE_models <- results$CATE_models              # CATE model used in each fold

```

Printing results will display estimates and 95% confidence intervals for each AIPTW estimate and the overall treatment effect under the optimal treatment rule. It will also provide the proportion of the dataset treated under the OTR:

```
                               Results Aggregated Across k =  2  folds 
---------------------------------------------------------------------------------------------------- 
                         Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
---------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]       2.5607              0.0428              2.4768              2.6445              
E[Y(0) | d(Z) = 1]       1.7601              0.0369              1.6878              1.8324              
E[Y(0) - Y(1)]           0.4                 0.0227              0.3554              0.4445              

Proportion treated under OTR:  0.5012

```

Alternatively, nuisance models could be pre-fit for a given set of covariates `W`. This is helpful for cycling through multiple potential decision rules (multiple sets of `Z`).

```R
  
  # Nuisance model SuperLearner libraries
  sl.library.outcome <- c("SL.mean", "SL.glm", "SL.glm.interaction")      # libraries to use for outcome model
  sl.library.treatment <- c("SL.mean", "SL.glm", "SL.glm.interaction")    # libraries to use for treatment model
  sl.library.missingness <- c("SL.mean", "SL.glm", "SL.glm.interaction")  # libraries to use for missingness model
  
  # CATE model SuperLearner libraries
  sl.library.CATE <- c("SL.mean", "SL.glm", "SL.glm.interaction")
  
  # List of covariates to use to estimate Nuisance models
  W_list <- c("W1", "W2")
  
  # List of covariates to use to estimate CATE model 
  Z_list <- c("W1")
  
  decision_threshold <- 0  # clinically relevance threshold for treatment effect (>=0 if desired outcome Y, negative in undesirable)
  
  nuisance_output <- learn_nuisance(df = df,
                                    Y_name = "Y",
                                    A_name = "A",
                                    W_list = W_list,
                                    sl.library.outcome = sl.library.outcome,
                                    sl.library.treatment = sl.library.treatment,
                                    sl.library.missingness = sl.library.missingness,
                                    outcome_type = "gaussian",
                                    k_folds = 2,
                                    ps_trunc_level = 0.01)
  
  nuisance_models <- nuisance_output$nuisance_models
  k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
  
  Z_lists <- list(c("W1"), c("W2"), c("W1", "W2"))
  results_list <- vector(mode = "list", length = length(Z_lists))
  
  for(i in 1:length(Z_lists)){
    Z_list <- Z_lists[[i]]
    
    results <- estimate_OTR(df = df,
                          Y_name = "Y",
                          A_name = "A",
                          W_list = W_list,
                          Z_list = Z_list,
                          sl.library.CATE = sl.library.CATE,
                          nuisance_models = nuisance_models,
                          k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                          threshold = 0,
                          k_folds = 2,
                          ps_trunc_level = 0.01,
                          outcome_type = "gaussian")
                          
    results_list[[i]] <- results
    
  }
  
```
