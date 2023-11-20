
# `drotr`

 MODIFY FOR NEW REPO/PACKAGE
## Scripts:
**1. estimate_EYd.R** - This is the primary R script for testing treatment rules in the ABCD study.

**2. simulation_EYd.R** - This script is used to run simulations during testing.

**3. generate_sim_data.R** - This script is used to generate simulated data. These simulated datasets were used during testing.

## How to use:

### `estimate_EYd.R`

The `estimate_EYd.R` script contains functions necessary to estimate the expected outcome for participants in the ABCD study. The main function you'll be working with is `estimate_EYd`.

**Usage**:

Example with simulation data

```R
  # Example usage of estimate_EYd function
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
  
  output_est <- estimate_EYd(df = df,
                             Y_name = "Y",
                             A_name = "A",
                             Z_list = Z_list, 
                             sl.library.outcome = sl.library, 
                             sl.library.treatment = sl.library, 
                             sl.library.missingness = sl.library,
                             sl.library.CATE = sl.library, 
                             threshold = decision_threshold,
                             k_folds = 2,
                             outcome_type = "gaussian", 
                             pwd = "~/Documents/Research/abcd_cate")
  
    overall_results <- output_est[[1]]                  # Overall results dataframe
    results_by_fold_aiptw_1 <- output_est[[2]]          # Results for each of k folds AIPTW, a = 1 
    results_by_fold_aiptw_0 <- output_est[[3]]          # Results for each of k folds AIPTW, a = 0
    results_by_fold_treatment_effect <- output_est[[4]] # Results for each of k folds treatment effect
    cate_models_each_fold <- output_est[[5]]            # CATE model used in each fold
    decision_each_participant <- output_est[[6]]        # Simulation data with decision made for each participant
    
    prop_treated <- mean(decision_each_participant$decision) #Proportion treated by rule

   ```
   
   Example with ABCD Data
   
   ```R
   # Example usage of estimate_EYd function
   # Treatment rule: Single pathogen quantity (shigella)
   
   # Libraries to try in each step
   sl.library.outcome <- c("SL.glm", "SL.ranger", "SL.lasso", "SL.earth", "SL.ManuallyDefinedLogRegression")
   sl.library.treatment <- c("SL.mean", "SL.ManuallyDefinedModels")
   sl.library.missingness <- c("SL.mean", "SL.ManuallyDefinedModels")
   sl.library.CATE <- c("SL.ManuallyDefinedModels")
   
   output_est <- estimate_EYd(Y_name = "day3diar",
                              A_name = "an_grp_01",
                              Z_list = "shigella_new", 
                              sl.library.outcome = sl.library.outcome, 
                              sl.library.treatment = sl.library.treatment, 
                              sl.library.missingness = sl.library.missingness,
                              sl.library.CATE = sl.library.CATE, 
                              threshold = 0.05,
                              k_folds = 2,
                              outcome_type = "binomial", 
                              pwd = "~/Documents/Research/abcd_cate")
                             
    overall_results <- output_est[[1]]                  # Overall results dataframe
    results_by_fold_aiptw_1 <- output_est[[2]]          # Results for each of k folds AIPTW, a = 1 
    results_by_fold_aiptw_0 <- output_est[[3]]          # Results for each of k folds AIPTW, a = 0
    results_by_fold_treatment_effect <- output_est[[4]] # Results for each of k folds treatment effect
    cate_models_each_fold <- output_est[[5]]            # CATE model used in each fold
    decision_each_participant <- output_est[[6]]        # Simulation data with decision made for each participant
    
    prop_treated <- mean(decision_each_participant$decision) #Proportion treated by rule
    
   ```
   
