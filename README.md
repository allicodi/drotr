
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
  
```

The function `estimate_OTR` will assign treatment to all observations in `df` with a clinically relevant individual level treatment effect. It will also estimate the treatment effect of those assigned treatment by the OTR. 

```R
  
  # Nuisance model SuperLearner libraries
  sl.library.outcome <- c("SL.glm", "SL.glm.interaction")      # libraries to use for outcome model
  sl.library.treatment <- c("SL.glm", "SL.glm.interaction")    # libraries to use for treatment model
  sl.library.missingness <- c("SL.mean", "SL.glm", "SL.glm.interaction")  # libraries to use for missingness model
  
  # CATE model SuperLearner libraries
  sl.library.CATE <- c("SL.mean", "SL.glm", "SL.glm.interaction")
  
  # List of covariates to use to estimate Nuisance models
  W_list <- c("W1", "W2")
  
  # List of covariates to use to estimate CATE model 
  Z_list <- c("W1")
  
  # Clinically relevant threshold for treatment effect
  # Zero will default to positive (desirable) outcome variable Y. If Y is undesirable outcome, set to "-0"
  decision_threshold <- "0"  
  
  otr_estimate <- estimate_OTR(df = df,
                                Y_name = "Y",
                                A_name = "A",
                                W_list = W_list,
                                Z_list = Z_list,
                                id_name = NULL, 
                                sl.library.CATE = sl.library.CATE,
                                sl.library.outcome = sl.library.outcome,
                                sl.library.treatment = sl.library.treatment,
                                sl.library.missingness = sl.library.missingness,
                                threshold = decision_threshold,
                                k_folds = 2,
                                ps_trunc_level = 0.01,
                                outcome_type = "gaussian")
    
  otr_estimate

```

Printing the `otr_estimate` results will display estimates and 95% confidence intervals for each AIPTW estimate, the proportion of patients treated under the decision rule, the subgroup treatment effect, and the overall treatment effect:

```
                               Results for  threshold =  0  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]            2.5842              0.0438              2.4984              2.67                
E[Y(0) | d(Z) = 1]            1.7756              0.0369              1.7032              1.8479              
E[d(Z) = 1]                   0.519               0.0071              0.5051              0.5329              
E[Y(d) - Y(0) | d(Z) = 1]     0.8087              0.043               0.7243              0.893               
E[Y(d) - Y(0)]                0.4195              0.023               0.3743              0.4646              

Covariates used in decision rule:  W1

```

A user could input a vector of different thresholds for which to identify treatment effects. 

A positive threshold indicates a desirable outcome (treat if CATE > threshold), and a negative threshold indicates an undesirable outcome (treat if CATE< threshold). A threshold of 0 defaults to positive outcome, but "-0" (as a string) can be specified to treat if CATE < threshold. 

```
  # Thresholds to test for treatment effect
  # Zero will default to positive (desirable) outcome variable Y. If Y is undesirable outcome, set to "-0"
  # Let's now assume Y is undesirable outcome in our simulation data
  decision_thresholds <- c("-0", "-0.05", "-0.10")  
  
  otr_estimate <- estimate_OTR(df = df,
                                Y_name = "Y",
                                A_name = "A",
                                W_list = W_list,
                                Z_list = Z_list,
                                id_name = NULL, 
                                sl.library.CATE = sl.library.CATE,
                                sl.library.outcome = sl.library.outcome,
                                sl.library.treatment = sl.library.treatment,
                                sl.library.missingness = sl.library.missingness,
                                threshold = decision_thresholds,
                                k_folds = 2,
                                ps_trunc_level = 0.01,
                                outcome_type = "gaussian")
    
  otr_estimate

```
```
                      Results for  threshold =  -0  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]            2.5863              0.0439              2.5004              2.6723              
E[Y(0) | d(Z) = 1]            1.7706              0.0368              1.6985              1.8426              
E[d(Z) = 1]                   0.5182              0.0071              0.5043              0.5321              
E[Y(d) - Y(0) | d(Z) = 1]     0.8158              0.0429              0.7317              0.8999              
E[Y(d) - Y(0)]                0.4227              0.023               0.3776              0.4678              

Covariates used in decision rule:  W1 

                      Results for  threshold =  -0.05  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]            2.5409              0.0431              2.4565              2.6253              
E[Y(0) | d(Z) = 1]            1.7522              0.0362              1.6812              1.8232              
E[d(Z) = 1]                   0.5348              0.0071              0.521               0.5486              
E[Y(d) - Y(0) | d(Z) = 1]     0.7887              0.0422              0.706               0.8714              
E[Y(d) - Y(0)]                0.4218              0.0233              0.3762              0.4673              

Covariates used in decision rule:  W1 

                      Results for  threshold =  -0.10  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]            2.477               0.0425              2.3937              2.5604              
E[Y(0) | d(Z) = 1]            1.7291              0.0356              1.6592              1.799               
E[d(Z) = 1]                   0.5552              0.007               0.5414              0.569               
E[Y(d) - Y(0) | d(Z) = 1]     0.748               0.0415              0.6665              0.8294              
E[Y(d) - Y(0)]                0.4152              0.0237              0.3688              0.4616              

Covariates used in decision rule:  W1 
```

Alternatively, nuisance models could be pre-fit for a given set of covariates `W`. This is helpful for cycling through multiple potential decision rules (multiple sets of `Z`).

```R
  
  # Nuisance model SuperLearner libraries
  sl.library.outcome <- c("SL.glm", "SL.glm.interaction")      # libraries to use for outcome model
  sl.library.treatment <- c("SL.glm", "SL.glm.interaction")    # libraries to use for treatment model
  sl.library.missingness <- c("SL.mean", "SL.glm", "SL.glm.interaction")  # libraries to use for missingness model
  
  # CATE model SuperLearner libraries
  sl.library.CATE <- c("SL.glm", "SL.glm.interaction")
  
  # List of covariates to use to estimate Nuisance models
  W_list <- c("W1", "W2")
  
  # Clinically relevant threshold for treatment effect
  # Zero will default to positive (desirable) outcome variable Y. If Y is undesirable outcome, set to "-0"
  decision_threshold <- 0  
  
  nuisance_output <- learn_nuisance(df = df,
                                    Y_name = "Y",
                                    A_name = "A",
                                    W_list = W_list,
                                    id_name = NULL, 
                                    sl.library.outcome = sl.library.outcome,
                                    sl.library.treatment = sl.library.treatment,
                                    sl.library.missingness = sl.library.missingness,
                                    outcome_type = "gaussian",
                                    k_folds = 2,
                                    ps_trunc_level = 0.01)
  
  nuisance_models <- nuisance_output$nuisance_models
  k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
  
  # List of different sets of covariates to use to estimate CATE model 
  Z_lists <- list(c("W1"), c("W2"), c("W1", "W2"))
  results_list <- vector(mode = "list", length = length(Z_lists))
  
  for(i in 1:length(Z_lists)){
    Z_list <- Z_lists[[i]]
    
    results <- estimate_OTR(df = df,
                          Y_name = "Y",
                          A_name = "A",
                          W_list = W_list,
                          Z_list = Z_list,
                          id_name = NULL,
                          sl.library.CATE = sl.library.CATE,
                          nuisance_models = nuisance_models,
                          k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                          threshold = 0,
                          k_folds = 2,
                          ps_trunc_level = 0.01,
                          outcome_type = "gaussian")
                          
    results_list[[i]] <- results
    
  }
  
  results_list
  
```

```
[[1]]
                      Results for  threshold =  0  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]            2.5873              0.0438              2.5013              2.6732              
E[Y(0) | d(Z) = 1]            1.7814              0.0366              1.7097              1.8532              
E[d(Z) = 1]                   0.5186              0.0071              0.5048              0.5324              
E[Y(d) - Y(0) | d(Z) = 1]     0.8058              0.043               0.7215              0.8901              
E[Y(d) - Y(0)]                0.417               0.023               0.372               0.462               

Covariates used in decision rule:  W1 


[[2]]
                      Results for  threshold =  0  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]            1.5859              0.0493              1.4893              1.6826              
E[Y(0) | d(Z) = 1]            1.5175              0.0341              1.4507              1.5842              
E[d(Z) = 1]                   0.704               0.0049              0.6944              0.7136              
E[Y(d) - Y(0) | d(Z) = 1]     0.0685              0.0436              -0.017              0.154               
E[Y(d) - Y(0)]                0.0231              0.0276              -0.0309             0.0771              

Covariates used in decision rule:  W2 


[[3]]
                      Results for  threshold =  0  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
E[Y(d) | d(Z) = 1]            2.6106              0.0437              2.5249              2.6962              
E[Y(0) | d(Z) = 1]            1.8021              0.0366              1.7304              1.8738              
E[d(Z) = 1]                   0.5168              0.0071              0.503               0.5306              
E[Y(d) - Y(0) | d(Z) = 1]     0.8084              0.0431              0.7239              0.893               
E[Y(d) - Y(0)]                0.4168              0.0229              0.3718              0.4618              

Covariates used in decision rule:  W1, W2 
```

Aggregated and individual level decisions by threshold and fold can be accessed through the results object. Models for outcome, treatment, missingness, and CATE by fold are also returned in the results object. 


