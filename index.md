# R/`drotr`

> Doubly robust optimal treatment rule estimation

**Author:** Allison Codi

## Description

`drotr` is a package that uses doubly-robust methods to estimate optimal
treatment rules (OTRs).

We define OTRs through the use of conditional average treatment effects
(CATEs). The CATE can be interpreted as the expected difference in
outcome under treatment vs placebo. A doubly-robust learner is used to
estimate the CATE. All individuals with CATE estimates that exceed a
specified threshold ‘t’ are treated under the OTR.

Once observations are assigned treatment under the OTR, we can use
various estimands to describe the rule: (1) the proportion treated, (2)
the overall Average Treatment effect under the Rule ($`ATR_d`$), (3) the
Average Treatment effect in the subgroup Recommended Treatment under the
rule ($`ATRT_d`$), (4) the Average Treatment effect in the subgroup Not
Recommended Treatment under the rule ($`ATNRT_d`$), and (5) the
difference between treatment recommendation subgroups
($`ATRT_d - ATNRT_d`$). These are estimated using Augmented Inverse
Probability of Treatment Weight (AIPTW) estimators.

The full procedure is implemented using nested cross-validation and
Super Learning.

## Installation:

A developmental release may be installed from GitHub via devtools with:

`devtools::install_github("allicodi/drotr")`

## Usage:

Here we demonstrate calls to `drotr` using a simulated data set of
length-for-age z-score 90 days post enrollment for a diarrheal
antibiotics trial.

``` r


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

# List of libraries to use for nuisance and CATE models
sl.library.outcome <- c("SL.glm", "SL.ranger", "SL.earth")
sl.library.treatment <- c("SL.mean", "SL.glm")
sl.library.missingness <- c("SL.mean", "SL.glm")
sl.library.CATE <- c("SL.glm", "SL.ranger", "SL.earth")

set.seed(12345)

results_host <- estimate_OTR(df = abcd_data,
                             Y_name = Y, 
                             A_name = A, 
                             W_list = W, 
                             Z_list = Z, 
                             id_name = "pid",
                             sl.library.outcome = sl.library.outcome,
                             sl.library.treatment = sl.library.treatment,
                             sl.library.missingness = sl.library.missingness,
                             sl.library.CATE = sl.library.CATE,
                             threshold = t,
                             k_folds = 10,
                             outcome_type = "gaussian")
                             
print(results_host)
                             
```

                                   Results for  threshold =  0.05  Aggregated Across k =  10  folds
    ---------------------------------------------------------------------------------------------------------------------------------------
                                                                Estimate            Standard Error      95% CI: Lower       95% CI: Upper
    ---------------------------------------------------------------------------------------------------------------------------------------
    E[Y(1) | d(Z) = 1]                                          -1.659              0.0266              -1.7112             -1.6069
    E[Y(0) | d(Z) = 1]                                          -1.7054             0.026               -1.7564             -1.6543
    E[d(Z) = 1]                                                 0.4507              0.0061              0.4388              0.4626
    E[Y(1) - Y(0) | d(Z) = 1]                                   0.0463              0.0209              0.0053              0.0874
    E[Y(1) - Y(0) | d(Z) = 0]                                   0.0337              0.0193              -0.0041             0.0716
    E[Y(d) - Y(0)]                                              0.0191              0.0094              6e-04               0.0376
    E[Y(1) - Y(0) | d(Z) = 1] - E[Y(1) - Y(0) | d(Z) = 0]       0.0126              0.0285              -0.0432             0.0684
    E[Y(d)]                                                     -1.6588             0.0172              -1.6926             -1.625
    E[Y]                                                        -1.6587             0.0156              -1.6892             -1.6282
    E[Y(d) - Y]                                                 -1e-04              0.0071              -0.0139             0.0138

    Covariates used in decision rule:  avemuac, wfazscore, wflzscore, lfazscore, dy1_ant_sex, agemchild, an_ses_quintile 

Additional information regarding rule interpretation, comparison, and
averaging across seeds can be found in `vignettes/drotr`

## Issues

If you encounter any bugs or have feature requests, please [file an
issue](https://github.com/allicodi/drotr/issues).

## Citation

After using the `drotr` R package, please cite the following:

    @Manual{drotr_package,
      title = {drotr: Doubly-robust optimal treatment rule estimation},
      author = {Allison Codi},
      note = {R package version 1.1.0}
    }
