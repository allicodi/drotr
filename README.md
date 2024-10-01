
# `drotr`

**Author:** Allison Codi

## Description

`drotr` is a package that uses doubly-robust methods to estimate optimal treatment rules (OTRs).

We define OTRs through the use of conditional average treatment effects (CATEs). The CATE can be interpreted as the expected difference in outcome under treatment vs placebo. A doubly-robust learner is used to estimate the CATE. All individuals with CATE estimates that exceed a specified threshold 't' are treated under the OTR. 

Once observations are assigned treatment under the OTR, we can use various estimands to describe the rule: (1) the proportion treated, (2) the overall Average Treatment effect under the Rule (ATR_d), (3) the Average Treatment effect in the subgroup Recommended Treatment under the rule (ATRT_d), (4) the Average Treatment effect in the subgroup Not Recommended Treatment under the rule (ATNRT_d), and (5) the difference between treatment recommendation subgroups (ATRT_d - ATNRT_d).
These are estimated using Augmented Inverse Probability of Treatment Weight (AIPTW) estimators.

The full procedure is implemented using nested cross-validation and Super Learning.

## Installation:

A developmental release may be installed from GitHub via devtools with:

```devtools::install_github("allicodi/drotr")```

## Usage:

Here we demonstrate calls to `drotr` using a simulated data set of length-for-age z-score 90 days post enrollment for a diarrheal antibiotics trial.

```R

library(drotr)

# 1. Load data
data(abcd_data)

# 2. Set input parameters

# outcome variable
Y <- "lazd90"

# treatment variable
A <- "an_grp_01"

# covariates 
W <- c("rotavirus_new", "rotavirus_bin", "norovirus_new", "norovirus_bin", "adenovirus_new",
      "adenovirus_bin", "sapovirus_new","sapovirus_bin", "astrovirus_new", "astrovirus_bin",
      "st_etec_new", "st_etec_bin", "shigella_new", "shigella_bin", "campylobacter_new",
      "campylobacter_bin", "tepec_new", "tepec_bin", "v_cholerae_new", "v_cholerae_bin",
      "salmonella_new", "salmonella_bin", "cryptosporidium_new", "cryptosporidium_bin",
      "dy1_scrn_vomitall", "dy1_scrn_lstools", "dy1_scrn_sstools", "dy1_scrn_diardays",
      "dy1_scrn_dehydr", "avemuac", "wfazscore", "lfazscore", "wflzscore", "site",
      "dy1_ant_sex", "agemchild", "an_ses_quintile", "an_tothhlt5", "month_en", "rotaseason")

# subset of covariates to use to make rule
Z <- c("avemuac", "wfazscore", "wflzscore", "lfazscore", "dy1_ant_sex", "agemchild", "an_ses_quintile")

# treatment threshold
t <- 0.05
```

## Issues

If you encounter any bugs or have feature requests, please [file an issue](https://github.com/allicodi/drotr/issues).


## Citation

After using the `drotr` R package, please cite the following:

@Manual{drotr_package,
  title = {drotr: Doubly-robust optimal treatment rule estimation},
  author = {Allison Codi},
  note = {R packag version 1.0.0}
}


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
                                k_folds = 5,
                                ps_trunc_level = 0.01,
                                outcome_type = "gaussian")
    
  otr_estimate

```

Printing the `otr_estimate` results will display estimates and 95% confidence intervals for the following estimands:

1. E[Y(1) | d(Z) = 1] - estimated outcome if treated among those recommended treatment under optimal treatment rule
2. E[Y(0) | d(Z) = 1] - estimated outcome if NOT treated among those recommended treatment under optimal treatment rule
3. E[d(Z) = 1] - estimated proportion of patients treated under optimal treatment rule
4. E[Y(1) - Y(0) | d(Z) = 1] - estimated treatment effect in subgroup of patients recommended treatment under optimal treatment rule
5. E[Y(1) - Y(0) | d(Z) = 0] - estimated treatment effect in subgroup of patients NOT recommended treatment under optimal treatment rule
6. E[Y(d) - Y(0) ] - overall treatment effect of the optimal treatment rule
7. E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0] - comparison of treatment effects between subgroups of those recommended and not recommended treatment

```
                                Results for  threshold =  0  Aggregated Across k =  5  folds 
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          2.5982              0.043               2.5139              2.6825              
E[Y(0) | d(Z) = 1]                                          1.8235              0.0379              1.7491              1.8978              
E[d(Z) = 1]                                                 0.4822              0.0071              0.4683              0.4961              
E[Y(1) - Y(0) | d(Z) = 1]                                   0.7747              0.0436              0.6893              0.8601              
E[Y(1) - Y(0) | d(Z) = 0]                                   -0.7899             0.0442              -0.8764             -0.7033             
E[Y(d) - Y(0)]                                              0.3734              0.0217              0.3308              0.416               
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       1.5645              0.062               1.4429              1.6861              

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
                                k_folds = 5,
                                ps_trunc_level = 0.01,
                                outcome_type = "gaussian")
    
  otr_estimate

```
```
                               Results for  threshold =  -0  Aggregated Across k =  5  folds 
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          -0.5736             0.0421              -0.6561             -0.491              
E[Y(0) | d(Z) = 1]                                          0.2306              0.0378              0.1566              0.3047              
E[d(Z) = 1]                                                 0.5168              0.0071              0.503               0.5306              
E[Y(1) - Y(0) | d(Z) = 1]                                   -0.8042             0.0441              -0.8905             -0.7178             
E[Y(1) - Y(0) | d(Z) = 0]                                   0.7833              0.0435              0.6979              0.8686              
E[Y(d) - Y(0)]                                              -0.4148             0.0235              -0.4609             -0.3687             
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       -1.5875             0.0619              -1.7089             -1.4661             

Covariates used in decision rule:  W1 

                               Results for  threshold =  -0.05  Aggregated Across k =  5  folds 
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          -0.6187             0.0427              -0.7024             -0.5349             
E[Y(0) | d(Z) = 1]                                          0.2001              0.0384              0.1249              0.2753              
E[d(Z) = 1]                                                 0.499               0.0071              0.4851              0.5129              
E[Y(1) - Y(0) | d(Z) = 1]                                   -0.8187             0.0449              -0.9068             -0.7307             
E[Y(1) - Y(0) | d(Z) = 0]                                   0.7408              0.0429              0.6566              0.825               
E[Y(d) - Y(0)]                                              -0.4077             0.0232              -0.4531             -0.3623             
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       -1.5595             0.0621              -1.6813             -1.4377             

Covariates used in decision rule:  W1 

                               Results for  threshold =  -0.10  Aggregated Across k =  5  folds 
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          -0.6813             0.0435              -0.7665             -0.5961             
E[Y(0) | d(Z) = 1]                                          0.1827              0.0393              0.1056              0.2598              
E[d(Z) = 1]                                                 0.4776              0.0071              0.4638              0.4914              
E[Y(1) - Y(0) | d(Z) = 1]                                   -0.864              0.0457              -0.9537             -0.7744             
E[Y(1) - Y(0) | d(Z) = 0]                                   0.7186              0.0421              0.636               0.8012              
E[Y(d) - Y(0)]                                              -0.4119             0.0227              -0.4564             -0.3673             
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       -1.5826             0.0622              -1.7045             -1.4608             

Covariates used in decision rule:  W1 
```

Alternatively, nuisance models could be pre-fit for a given set of covariates `W`. The `validRows` argument must be used to pass fold assignments through from step to step. 

Pre-fitting nuisance models may be helpful when cycling through multiple potential decision rules (multiple sets of `Z`).

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
  validRows <- nuisance_output$validRows
  
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
                          validRows = validRows,
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
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          2.6006              0.043               2.5162              2.6849              
E[Y(0) | d(Z) = 1]                                          1.8179              0.038               1.7434              1.8925              
E[d(Z) = 1]                                                 0.4838              0.0071              0.4699              0.4977              
E[Y(1) - Y(0) | d(Z) = 1]                                   0.7827              0.0437              0.697               0.8683              
E[Y(1) - Y(0) | d(Z) = 0]                                   -0.7946             0.0443              -0.8814             -0.7079             
E[Y(d) - Y(0)]                                              0.3787              0.0219              0.3358              0.4216              
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       1.5773              0.0622              1.4554              1.6992              

Covariates used in decision rule:  W1 


[[2]]
                               Results for  threshold =  0  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          -9e-04              0.0617              -0.1218             0.1199              
E[Y(0) | d(Z) = 1]                                          0.0366              0.048               -0.0574             0.1306              
E[d(Z) = 1]                                                 0.4226              0.0062              0.4104              0.4348              
E[Y(1) - Y(0) | d(Z) = 1]                                   -0.0376             0.0569              -0.1491             0.0739              
E[Y(1) - Y(0) | d(Z) = 0]                                   -0.0422             0.047               -0.1342             0.0499              
E[Y(d) - Y(0)]                                              -0.0122             0.0212              -0.0538             0.0294              
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       0.0046              0.0738              -0.14               0.1492              

Covariates used in decision rule:  W2 


[[3]]
                               Results for  threshold =  0  Aggregated Across k =  2  folds 
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          2.5928              0.0432              2.5081              2.6775              
E[Y(0) | d(Z) = 1]                                          1.8004              0.0382              1.7255              1.8754              
E[d(Z) = 1]                                                 0.48                0.0071              0.4661              0.4939              
E[Y(1) - Y(0) | d(Z) = 1]                                   0.7923              0.0438              0.7065              0.8782              
E[Y(1) - Y(0) | d(Z) = 0]                                   -0.7924             0.0441              -0.8789             -0.7059             
E[Y(d) - Y(0)]                                              0.3804              0.0218              0.3377              0.4231              
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       1.5848              0.0621              1.463               1.7066              

Covariates used in decision rule:  W1, W2 
```

The full results object (class `full_otr_results`) has the following structure:

```
results_object
|_results (class `otr_results`)  
  |_threshold = t1 `Results` object  
    |_aggregated_results  
    |_k_fold_results  
    |_decision_df
  |_ threshold = t2 `Results` object  
  |_ ...  
  |_ threshold = tn `Results` object  
  |_Z_list
|_nuisance_models  
  |_fold 1 `Nuisance` object  
    |_outcome_model  
    |_treatment_model  
    |_missingess_model  
  |_fold 2 `Nuisance` object  
  |_ ...  
  |_fold k `Nuisance` object  
|_CATE_models  
  |_ fold 1 CATE model  
  |_ fold 2 CATE model  
  |_ ... 
  |_ fold k CATE model  
|_Z_list  
```

Aggregated and individual level decisions by threshold and fold can be accessed through the results object. Models for outcome, treatment, missingness, and CATE by fold are also returned in the results object. More details about components of results object can be found in roxygen help documentation.

The compare.otr_results function can be used to compare treatment effects and subgroup effects across treatment rules and thresholds.

```
compare.otr_results(results_list[[1]], results_list[[2]], 0, "te", "te")
```

```
 Treatment Effect E[Y(d) - Y(0)]  for rule 1 at threshold =  0 
 vs 
 Treatment Effect E[Y(d) - Y(0)]  for rule 2 at threshold =  0 
--------------------------------------------------------------------------------------------------------- 
                              Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------- 
Rule 1 - Rule 2               0.3908              0.0221              0.3474              0.4342              

 Rule 1: Z =  W1
 Rule 2: Z =  W2
 
```

It may be helpful to repeat analyses across multiple seeds to examine average performance. 
The `average_across_seeds` helper function can be used to assist with averaging results. 

```

results_list <- vector(mode = "list", length = 5)

for(seed in 1:5){

  set.seed(seed)
  results_list[[seed]] <- estimate_OTR(df = df,
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
                                k_folds = 5,
                                ps_trunc_level = 0.01,
                                outcome_type = "gaussian")

}

average_across_seeds(results_list, -0.05)

```

```
                                          Average results across n =  5  seeds for threshold  -0.05 
--------------------------------------------------------------------------------------------------------------------------------------- 
                                                            Estimate            Standard Error      95% CI: Lower       95% CI: Upper       
--------------------------------------------------------------------------------------------------------------------------------------- 
E[Y(1) | d(Z) = 1]                                          -0.6031             0.0431              -0.6875             -0.5187             
E[Y(0) | d(Z) = 1]                                          0.2253              0.0386              0.1498              0.3009              
E[d(Z) = 1]                                                 0.4978              0.0071              0.4839              0.5117              
E[Y(1) - Y(0) | d(Z) = 1]                                   -0.8285             0.0451              -0.9168             -0.7402             
E[Y(1) - Y(0) | d(Z) = 0]                                   0.7448              0.0427              0.661               0.8286              
E[Y(d) - Y(0)]                                              -0.4123             0.0232              -0.4578             -0.3668             
E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       -1.5733             0.0621              -1.695              -1.4515             

Covariates used in decision rule:  W1, W2 
```
