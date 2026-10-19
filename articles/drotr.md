# drotr: Learning and evaluating threshold-based treatment rules

Abstract

The drotr package implements a cross-validation-based framework for
learning, evaluating, and comparing individualized treatment rules.
Rules are defined by thresholding estimates of the conditional average
treatment effect (CATE) to recommend treatment only when the expected
benefit reaches a prespecified clinically meaningful magnitude. The
framework allows rules based on different sets of baseline covariates to
be compared, facilitating evaluation of the potential value of
additional measurements such as diagnostic biomarkers. Nested
cross-validation, a doubly robust learner, and Super Learner are used to
develop the rules, with augmented inverse probability weighted
estimators used to evaluate outcomes and effects under each rule. This
vignette introduces the estimands and demonstrates rule learning,
evaluation on held-out data, and comparison of treatment rules using
data simulated to resemble the AntiBiotics for Children with severe
Diarrhea (ABCD) trial.

## Introduction

Individualized treatment rules tailor recommendations to patient
characteristics. Such rules are useful when the benefit of treatment
varies across individuals, but treating everyone with any positive
expected benefit may be undesirable. For example, antibiotics can
benefit children with bacterial diarrheal illness, yet broad antibiotic
use can contribute to antimicrobial resistance. Treatment guidance
therefore needs to weigh expected clinical benefit against the
consequences of unnecessary treatment.

The `drotr` package implements the framework described in *Machine
learning for estimating and comparing clinical rules for treating
diarrheal illness with antibiotics* (Codi et al.). The approach focuses
on two practical questions: **who is expected to benefit enough to
justify treatment**, and **which patient characteristics are needed to
make useful treatment recommendations**? We define rules by thresholding
the conditional average treatment effect (CATE), allowing a prespecified
minimum magnitude of expected clinical benefit. We then evaluate and
compare rules based on different subsets of covariates, including
readily available clinical characteristics and potentially more costly
diagnostic measurements.

Treatment rules are learned using a doubly robust (DR) learner (Kennedy
2023), with Super Learner used for nuisance and CATE regression
estimation (Luedtke and Laan 2016; Montoya et al. 2023). An outer
cross-validation procedure separates rule development from effect
evaluation. On held-out data, augmented inverse probability of treatment
weighted (AIPTW), or one-step, estimators are used to estimate
subgroup-specific treatment effects, the proportion recommended
treatment, and population-level outcomes under each rule. The resulting
inference pertains to the *data-adaptive rules learned in the training
folds*, not to an independently specified fixed rule.

## Methods

### Notation and treatment rules

Let $`W \in \mathcal{W}`$ denote measured baseline covariates and let
$`Z \subseteq W`$ represent the subset of covariates used to make
treatment recommendations. Let $`A \in \{0,1\}`$ indicate binary
treatment, $`Y`$ the outcome of interest, and $`\Delta \in \{0,1\}`$ an
indicator that $`Y`$ is observed ($`\Delta=1`$). The observed data
consist of independent observations
$`O_i=(W_i,A_i,\Delta_i,\Delta_iY_i)`$, $`i=1,\ldots,n`$. For
$`a \in \{0,1\}`$, let $`Y(a,1)`$ denote the outcome that would be
observed under treatment $`a`$ without outcome missingness. We write
$`Y(a)`$ for $`Y(a,1)`$ below to match the package’s printed output.

For a particular set of input covariates $`Z`$, the CATE is

``` math
\psi(z)=E[Y(1)-Y(0)\mid Z=z].
```

A treatment rule $`d:\mathcal{Z}\rightarrow\{0,1\}`$ assigns treatment
based on $`Z`$. Rather than requiring treatment for everyone with any
expected benefit, we use a clinically meaningful threshold $`t`$. For a
**desirable outcome**, such as length-for-age z-score, the package
defines a rule that recommends treatment for an estimated CATE greater
than $`t>0`$:

``` math
d(Z)=I\{\psi(Z)>t\}.
```

For an **undesirable outcome**, such as the presence of diarrhea,
improvement is a reduction in the outcome; the rule recommends treatment
when the CATE is less than $`t<0`$:

``` math
d(Z)=I\{\psi(Z)<t\}.
```

The rule is determined by both the selected covariates $`Z`$ and the
threshold $`t`$. In
[`estimate_OTR()`](https://allicodi.github.io/drotr/reference/estimate_OTR.md),
positive thresholds specify desirable outcomes and negative thresholds
specify undesirable outcomes. For a zero threshold, provide `"+0"` or
`"-0"` to indicate the desired direction. The package uses strict
inequalities when applying the threshold.

### Estimands for evaluating a rule

Let $`Y(d)=d(Z)Y(1)+\{1-d(Z)\}Y(0)`$ denote the potential outcome if
everyone received treatment according to rule $`d`$, and let
$`p_d=P\{d(Z)=1\}`$ be the proportion recommended treatment. The primary
estimands include:

1.  **Average treatment effect among those recommended treatment
    (ATRT)**: $`ATRT_d=E[Y(1)-Y(0)\mid d(Z)=1]`$.
2.  **Average treatment effect among those not recommended treatment
    (ATNRT)**: $`ATNRT_d=E[Y(1)-Y(0)\mid d(Z)=0]`$.
3.  **Proportion recommended treatment**: $`p_d=P\{d(Z)=1\}`$.
4.  **Average treatment effect under the rule (ATR)**:
    $`ATR_d=E[Y(d)-Y(0)]`$, comparing use of the rule against a policy
    of treating nobody.
5.  **Difference in subgroup treatment effects**: $`ATRT_d-ATNRT_d`$,
    evaluating how well the rule distinguishes groups with different
    expected treatment benefits.

In particular,

``` math
ATR_d=p_d\,ATRT_d,
```

because people not recommended treatment receive the same treatment
under rule $`d`$ and under the no-treatment policy.

The package also reports **outcome levels** in addition to treatment
effects:

- $`E[Y(1)\mid d(Z)=1]`$ and $`E[Y(0)\mid d(Z)=1]`$: expected outcomes
  under treatment and no treatment among those the rule recommends
  treating.
- $`E[Y(d)]`$: expected population outcome if everyone follows rule
  $`d`$.
- $`E[Y]`$: expected outcome under the observed treatment-assignment
  process.
- $`E[Y(d)-Y]`$: expected difference between outcomes under rule $`d`$
  and the observed treatment-assignment process.

The last three estimates are included in `drotr` result objects from
version 0.2.0 and later. When an older result object lacks the new
quantities, its print method displays `--`; a computed but undefined
estimate remains `NA`.

### Estimation and cross-validation

We use outer $`K`$-fold cross-validation to separate the development of
each treatment rule from its evaluation. Each fold has a training
sample, used to learn a rule, and an independent held-out validation
sample, used to estimate that rule’s performance. The training sample is
further divided into $`V`$ inner folds to construct a DR learner for the
CATE (the package uses three inner folds). This separation reduces
overfitting and supports inference on the average performance of the
cross-validated, data-adaptive rules.

Within an outer training sample, the DR learner uses three nuisance
regressions:

``` math
\begin{aligned}
\mu(a,w) &= E[Y\mid A=a,\Delta=1,W=w],\\
\pi(a\mid w) &= P(A=a\mid W=w),\\
\gamma(1\mid a,w) &= P(\Delta=1\mid A=a,W=w).
\end{aligned}
```

We fit these regressions using Super Learner. For observation $`i`$ in
an inner validation fold $`v`$ of outer training fold $`k`$, define the
pseudo-outcome

``` math
\widehat D_{k,v}(O_i)=
\frac{(2A_i-1)\Delta_i}
{\widehat\pi_{k,v}(A_i\mid W_i)\widehat\gamma_{k,v}(1\mid A_i,W_i)}
\{Y_i-\widehat\mu_{k,v}(A_i,W_i)\}
+\widehat\mu_{k,v}(1,W_i)-\widehat\mu_{k,v}(0,W_i).
```

Regressing the pseudo-outcome on $`Z`$ in each inner validation fold
yields an estimate of the CATE. We average these estimates across the
$`V`$ inner folds to define the rule learned from outer training fold
$`k`$, denoted $`\widehat d_k`$. For a desirable outcome, for example,
$`\widehat d_k(Z)=I\{\widehat\psi_k(Z)>t\}`$. Under the relevant
consistency and rate conditions, the DR construction allows consistent
CATE estimation if the outcome regression is consistently estimated or
if both the treatment and missingness regressions are consistently
estimated (Kennedy 2023).

To evaluate a learned rule on its outer validation sample, we
re-estimate the nuisance regressions using the full outer training
sample. For $`a,a'\in\{0,1\}`$, define

``` math
\theta_k(a\mid a')=
E\!\left[E(Y\mid A=a,\Delta=1,W)\mid\widehat d_k(Z)=a'\right],
\qquad
p_k(a')=P\{\widehat d_k(Z)=a'\}.
```

Writing $`\mathcal I_k`$ for the held-out observations,
$`n_k=|\mathcal I_k|`$, and $`\widehat p_k(a')`$ for the observed
proportion assigned recommendation $`a'`$ in that sample, the one-step
estimator can be expressed as

``` math
\widehat\theta_k(a\mid a')=\frac{1}{n_k}\sum_{i\in\mathcal I_k}
\frac{I\{\widehat d_k(Z_i)=a'\}}{\widehat p_k(a')}
\left[
\frac{I(A_i=a,\Delta_i=1)}
{\widehat\pi_k(a\mid W_i)\widehat\gamma_k(1\mid a,W_i)}
\{Y_i-\widehat\mu_k(a,W_i)\}
+\widehat\mu_k(a,W_i)
\right].
```

The subgroup effects are obtained from
$`\theta_k(1\mid1)-\theta_k(0\mid1)`$ and
$`\theta_k(1\mid0)-\theta_k(0\mid0)`$; the rule effect additionally uses
the proportion treated. Standard errors for these effects and
comparisons are based on influence-function calculations and the delta
method. The fold-specific estimates are then aggregated over outer
validation folds. Because the rules vary across training samples,
resulting confidence intervals should be interpreted for the
corresponding cross-validated data-adaptive targets, rather than for a
prespecified rule.

## Usage

First, load the package. **The model-fitting examples below are shown
but not executed when the vignette is built** (`eval = FALSE`) because
nested cross-validation with several Super Learner candidates is
computationally intensive. The data setup and learner definitions do run
during the build. To reproduce an analysis interactively, execute the
relevant model-fitting chunks and their dependent output chunks in
order.

First, load the package:

``` r

library(drotr)
```

The example uses `abcd_data`, a simulated dataset distributed with
`drotr` and designed to resemble data from the AntiBiotics for Children
with severe Diarrhea (ABCD) randomized trial. In ABCD, children with
acute watery diarrhea were randomized to azithromycin or placebo; the
simulated dataset contains $`n=6{,}692`$ observations and 40 baseline
covariates. This illustration uses simulated, not individual-level
trial, data.

``` r

data(abcd_data)
head(abcd_data[c("pid", "an_grp_01", "lazd90", "agemchild", "shigella_new")])
#>   pid an_grp_01      lazd90 agemchild shigella_new
#> 1   1         0 -0.24423622  5.801993     1.975475
#> 2   2         0 -3.24702621 14.681659     0.000000
#> 3   3         0 -0.05373623  7.558416     0.000000
#> 4   4         0 -2.17004791 13.743810     0.000000
#> 5   5         1  0.17351546  9.848812     0.000000
#> 6   6         1 -2.21432910  9.799398     0.000000
```

We consider azithromycin (`an_grp_01`, coded 1 for treatment and 0 for
placebo) and length-for-age z-score 90 days after enrollment (`lazd90`).
Approximately 4% of the simulated outcomes are missing. We first
construct a **host-characteristics rule** using a subset of baseline
anthropometric, demographic, and socioeconomic covariates (`Z`) to make
recommendations, while allowing the nuisance models to use a broader set
(`W`). With $`t=0.05`$, treatment is recommended to children whose
estimated CATE exceeds a gain of 0.05 length-for-age z-score units.

Higher length-for-age z-scores are desirable, so we use a positive
threshold. For an undesirable binary outcome such as diarrhea at day
three, a negative threshold would instead identify children whose
estimated probability of the outcome decreases sufficiently under
treatment. At zero, use the character string `"+0"` (desirable outcome)
or `"-0"` (undesirable outcome) to specify the direction.

### Inputs

Let’s set the parameters described above to be passed into the function:

``` r

Y <- "lazd90"

A <- "an_grp_01"

W <- c("rotavirus_new", "rotavirus_bin", "norovirus_new", "norovirus_bin", "adenovirus_new",
      "adenovirus_bin", "sapovirus_new","sapovirus_bin", "astrovirus_new", "astrovirus_bin",
      "st_etec_new", "st_etec_bin", "shigella_new", "shigella_bin", "campylobacter_new",
      "campylobacter_bin", "tepec_new", "tepec_bin", "v_cholerae_new", "v_cholerae_bin",
      "salmonella_new", "salmonella_bin", "cryptosporidium_new", "cryptosporidium_bin",
      "dy1_scrn_vomitall", "dy1_scrn_lstools", "dy1_scrn_sstools", "dy1_scrn_diardays",
      "dy1_scrn_dehydr", "avemuac", "wfazscore", "lfazscore", "wflzscore", "site",
      "dy1_ant_sex", "agemchild", "an_ses_quintile", "an_tothhlt5", "month_en", "rotaseason")

Z_host <- c("avemuac", "wfazscore", "wflzscore", "lfazscore", "dy1_ant_sex", "agemchild", "an_ses_quintile")

t <- 0.05
```

Next, we need to choose which libraries we are interested in using to
fit our nuisance models. The SuperLearner package comes with wrappers
for common prediction algorithms (see
[`SuperLearner::listWrappers()`](https://rdrr.io/pkg/SuperLearner/man/listWrappers.html)
for full list of available wrappers). We will choose several for our
outcome model that should be able to capture the complex form of the
outcome model.

``` r

sl.library.outcome <- c("SL.glm", "SL.ranger", "SL.earth")
```

Because treatment was randomized 1:1 in the motivating ABCD trial, a
constant-probability learner such as `SL.mean` is a natural choice for
the treatment propensity model in this example.

``` r

sl.library.treatment <- c("SL.mean")
```

We also model the probability that an outcome is observed, conditional
on treatment and measured covariates. The package supports custom Super
Learner wrappers; here, we illustrate a candidate regression for
missingness based on site and enrollment month, alongside a
constant-probability model. These two covariates are motivated by
patterns of outcome measurement in the ABCD trial.

``` r

SL.missing.1 <- function(Y, X, newX, family, ...){
  sl.missing.1_fit <- glm(Y ~ site + month_en,
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.missing.1_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.missing.1 = sl.missing.1_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.missing.1"
  # return the output
  return(out)
}
predict.SL.missing.1 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.missing.1, newdata = newdata, type="response")
  return(pred)
}
```

Include this wrapper among the candidate missingness models:

``` r

sl.library.missingness <- c("SL.mean", "SL.missing.1")
```

Lastly, specify candidate learners for regressing the DR pseudo-outcome
on the treatment-rule covariates $`Z`$:

``` r

sl.library.CATE <- c("SL.glm", "SL.ranger", "SL.earth")
```

### Estimate OTR

Call
[`estimate_OTR()`](https://allicodi.github.io/drotr/reference/estimate_OTR.md)
to fit the nuisance and CATE regressions, learn the threshold-based
rule, and estimate the corresponding effects. We identify participants
with `pid`, use ten outer cross-validation folds (`k_folds = 10`), and
specify a continuous (`"gaussian"`) outcome. This example can be
computationally intensive depending on the learners and available
hardware.

``` r

set.seed(12345)

results_host <- estimate_OTR(df = abcd_data,
                             Y_name = Y, 
                             A_name = A, 
                             W_list = W, 
                             Z_list = Z_host, 
                             id_name = "pid",
                             sl.library.outcome = sl.library.outcome,
                             sl.library.treatment = sl.library.treatment,
                             sl.library.missingness = sl.library.missingness,
                             sl.library.CATE = sl.library.CATE,
                             threshold = t,
                             k_folds = 10,
                             outcome_type = "gaussian")
```

As an alternative to the preceding **one-call** workflow, the analysis
can be split into two stages: (1) estimate and retain nuisance
regressions and fold assignments, then (2) estimate one or more
treatment rules using these fits. You need only run one of the two
workflows. Reusing nuisance fits is particularly helpful when comparing
rules using different sets of $`Z`$ or thresholds, because the rule
evaluations can share fold assignments. We use this staged workflow for
the comparison below. First, call
[`learn_nuisance()`](https://allicodi.github.io/drotr/reference/learn_nuisance.md):

``` r

set.seed(12345)
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = Y,
                                  A_name = A, 
                                  W_list = W,
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  k_folds = 10,
                                  outcome_type = "gaussian")
```

The returned object includes fold-specific nuisance model fits
(`nuisance_models`), outer fold assignments and DR pseudo-outcomes
(`k_fold_assign_and_CATE`), and inner validation indices (`validRows`).
Supply them to
[`estimate_OTR()`](https://allicodi.github.io/drotr/reference/estimate_OTR.md)
to avoid re-fitting the nuisance regressions.

``` r

nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

results_host <- estimate_OTR(df = abcd_data,
                             Y_name = Y, 
                             A_name = A, 
                             W_list = W, 
                             Z_list = Z_host, 
                             id_name = "pid",
                             sl.library.CATE = sl.library.CATE,
                             nuisance_models = nuisance_models,
                             k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                             validRows = validRows,
                             threshold = t,
                             k_folds = 10,
                             outcome_type = "gaussian")
```

Once the staged workflow has been run, we can print the estimates using
`print`:

``` r

print(results_host)
```

When run, the printed table reports estimates, standard errors, and 95%
confidence intervals. In this example, a **positive** effect indicates
improved length-for-age z-score under azithromycin. The rows have the
following interpretations:

- **$`E[Y(1)\mid d(Z)=1]`$ and $`E[Y(0)\mid d(Z)=1]`$:** expected day-90
  length-for-age z-scores if children recommended treatment were treated
  or untreated, respectively.
- **$`E[d(Z)=1]`$:** the proportion of children who would be recommended
  azithromycin under the host-characteristics rule.
- **$`E[Y(1)-Y(0)\mid d(Z)=1]`$ ($`ATRT_d`$):** treatment benefit in the
  subgroup selected for treatment. This measures the expected clinical
  benefit among those recommended treatment.
- **$`E[Y(1)-Y(0)\mid d(Z)=0]`$ ($`ATNRT_d`$):** treatment benefit in
  the subgroup not selected for treatment.
- **$`E[Y(d)-Y(0)]`$ ($`ATR_d`$):** expected population-level
  improvement from implementing the rule, compared with treating nobody.
- **$`ATRT_d-ATNRT_d`$:** difference in the average treatment benefit
  between the two recommendation subgroups. This evaluates whether the
  rule differentiates children by expected benefit.
- **$`E[Y(d)]`$:** expected population-average day-90 length-for-age
  z-score if treatment were assigned according to the rule.
- **$`E[Y]`$:** expected outcome under the observed treatment-assignment
  process.
- **$`E[Y(d)-Y]`$:** expected population-average change from the
  observed treatment-assignment process to implementing the rule.

One may wish to examine other output in the final `full_otr_results`
object. The `full_otr_results` object contains the following elements:

- **`results` (class `otr_results`)**  
  This component contains an object of class `Results` for each
  requested threshold.

  - **`threshold = t1` (Results object)**
    - `aggregated_results`: Estimated outcomes, effects, and standard
      errors aggregated across validation folds
    - `k_fold_results`: Fold-specific estimates and influence-function
      components
    - `decision_df`: Estimated CATEs and treatment recommendations in
      the held-out observations
    - `k_non_na`: Folds with estimable subgroup effects (excluding folds
      recommending treatment to everyone or nobody)
  - **`threshold = t2` (Results object)**
    - Same structure as `threshold = t1` if additional threshold(s)
      specified
  - **`...`**
    - Additional thresholds can be added in the same manner.

- **`nuisance_models`**  
  A list of nuisance models used in the analysis.

  - **`fold 1` (Nuisance object)**
    - `outcome_model`: Model for outcome.
    - `treatment_model`: Model for the treatment.
    - `missingness_model`: Model for missing data.
  - **`fold 2` (Nuisance object)**
    - Structure similar to `fold 1`.
  - **`...`**
    - Up to K nuisance models

- **`CATE_models`**  
  A list of CATE models generated during the analysis.

  - **`fold 1`**: CATE model results for fold 1.
  - **`fold 2`**: CATE model results for fold 2.
  - **`...`**: Up to k CATE models

- **`Z_list`**  
  List of covariates used for the CATE model

Model objects may be large, especially when using computationally
complex Super Learner libraries. To conserve storage while retaining
aggregated estimates and printing functionality, save the smaller
`otr_results` component (`results_host$results`) rather than the full
object if you no longer need the nuisance and CATE fits.

### Apply a learned OTR to held-out data

The
[`apply_OTR()`](https://allicodi.github.io/drotr/reference/apply_OTR.md)
function applies an already learned rule to another dataset. Rather than
re-learning the CATE, it obtains predictions from previously fitted CATE
models and applies the chosen threshold to assign treatment
recommendations. It then estimates the proportion recommended treatment,
subgroup treatment effects, and population-level outcomes under the
resulting rule using the outcome and treatment data from the evaluation
sample. This supports evaluating whether a rule developed in one dataset
is useful in another population or a previously reserved validation
sample.

The evaluation dataset must contain the **same CATE predictor
variables** (`Z_list`) as the development dataset, together with the
outcome (`Y_name`), treatment (`A_name`), and the covariates (`W_list`)
needed for outcome, treatment, and missingness estimation. These
variables should have compatible coding and measurement definitions
across datasets. For this evaluation to estimate causal effects, the
dataset needs observed treatment and outcome information (outcomes may
be missing for some observations). Predictor-only data can be used to
derive recommendations in principle, but do not suffice for the
effect-estimation workflow illustrated here.

To illustrate, we partition the simulated `abcd_data` into
nonoverlapping development and evaluation samples. We fit a rule in the
development sample using
[`estimate_OTR()`](https://allicodi.github.io/drotr/reference/estimate_OTR.md),
then supply its `CATE_models` to
[`apply_OTR()`](https://allicodi.github.io/drotr/reference/apply_OTR.md)
to evaluate the rule on the held-out sample. We use the same treatment
threshold, $`t=0.05`$, and the host-characteristic covariates defined
above (`Z_host`). Nuisance models are fitted within the evaluation data
using the supplied Super Learner libraries. This example is not executed
during vignette builds; to run it, execute its model-fitting and output
chunks.

``` r

set.seed(12345)

train_idx <- sample(seq_len(nrow(abcd_data)),
                          size = floor(0.8 * nrow(abcd_data)))
train_data <- abcd_data[train_idx, , drop = FALSE]
evaluation_data <- abcd_data[-train_idx, , drop = FALSE]

# Use Z_host as defined in the Inputs section; no earlier fitted model is needed.

# Learn the rule on the train sample only
rule_development <- estimate_OTR(
  df = train_data,
  Y_name = Y,
  A_name = A,
  W_list = W,
  Z_list = Z_host,
  id_name = "pid",
  sl.library.outcome = sl.library.outcome,
  sl.library.treatment = sl.library.treatment,
  sl.library.missingness = sl.library.missingness,
  sl.library.CATE = sl.library.CATE,
  threshold = t,
  k_folds = 5,
  outcome_type = "gaussian"
)

# Apply the learned CATE models to the held-out evaluation data
results_evaluation <- apply_OTR(
  df = evaluation_data,
  CATE_models = rule_development$CATE_models,
  Y_name = Y,
  A_name = A,
  W_list = W,
  Z_list = Z_host,
  id_name = "pid",
  threshold = t,
  sl.library.outcome = sl.library.outcome,
  sl.library.treatment = sl.library.treatment,
  sl.library.missingness = sl.library.missingness,
  k_folds = 5,
  outcome_type = "gaussian"
)

print(results_evaluation)
```

The returned `full_otr_results` object has the same general structure as
the output of
[`estimate_OTR()`](https://allicodi.github.io/drotr/reference/estimate_OTR.md).
For example, the treatment recommendations and model-specific CATE
predictions can be inspected through the `decision_df` component:

``` r

head(results_evaluation$results[[paste("threshold = ", t)]]$decision_df)
```

For the model-based option, the CATE predictions used to define the rule
are obtained by averaging predictions across the supplied trained CATE
models before applying the threshold. The reported estimates describe
the rule’s performance in the **evaluation sample**, rather than the
cross-validated performance of the rules learned within the training
sample. The interpretation of the estimands remains the same as above;
for example, $`E[Y(d)-Y(0)]`$ compares the rule with treating nobody,
whereas $`E[Y(d)-Y]`$ compares it with the observed treatment-assignment
process in the evaluation sample. The same causal identification and
estimation assumptions remain relevant.

Alternatively, if treatment decisions are already available as a binary
vector (one entry per observation in the evaluation data),
[`apply_OTR()`](https://allicodi.github.io/drotr/reference/apply_OTR.md)
can accept `d_pred` instead of `CATE_models`. In that case it evaluates
the supplied recommendations without predicting the CATE or applying the
threshold again. For evaluating a rule on genuinely external data,
nuisance model estimation must also be appropriate for the new
population; the example above fits these models on the evaluation
dataset rather than reusing models learned on the development data.

### Compare OTRs

An important goal of the framework is to compare rules built from
different input characteristics. For example, a rule using only readily
available host characteristics can be compared with a rule incorporating
pathogen diagnostics to investigate the potential value of measuring
additional biomarkers. Here we construct a second rule using *Shigella*
quantity and detection variables. Both rules use the same threshold,
outcome, nuisance fits, and cross-validation splits, so their comparison
uses aligned validation samples. The following code follows the staged
[`learn_nuisance()`](https://allicodi.github.io/drotr/reference/learn_nuisance.md)
example above; run those chunks first.

``` r

Z_shigella <- c("shigella_new", "shigella_bin")

results_shigella <- estimate_OTR(
  df = abcd_data,
  Y_name = Y,
  A_name = A,
  W_list = W,
  Z_list = Z_shigella,
  id_name = "pid",
  sl.library.CATE = sl.library.CATE,
  nuisance_models = nuisance_models,
  k_fold_assign_and_CATE = k_fold_assign_and_CATE,
  validRows = validRows,
  threshold = t,
  k_folds = 10,
  outcome_type = "gaussian"
)

print(results_shigella)
```

The
[`compare.otr_results()`](https://allicodi.github.io/drotr/reference/compare.otr_results.md)
function compares estimated effects between two rules while accounting
for their dependence from being evaluated on the same observations. Here
`"se"` specifies the subgroup effect $`ATRT_d`$ for each rule. Other
supported comparison targets include the population-level rule effect
(`"te"`) and the effect among those not recommended treatment
(`"se_dZ0"`).

``` r

compare.otr_results(res_rule1 = results_host,
                    res_rule2 = results_shigella,
                    threshold = t, 
                    rule1_comp = "se",
                    rule2_comp = "se")
```

For this call, the result is
$`ATRT_{d_{\mathrm{host}}}-ATRT_{d_{\mathrm{Shigella}}}`$. A negative
estimate indicates that the **average benefit among those recommended
treatment** is larger for the *Shigella*-based rule. Because the rules
may select different subgroups, this comparison is **not** a direct
contrast of average outcomes under the two policies, and it should not
be described as the difference in mean length-for-age z-scores between
the two selected groups. Its primary interpretation is the relative
ability of the two input sets to identify children with greater
treatment benefit. For population-level benefit, compare $`ATR_d`$ or
the expected outcomes $`E[Y(d)]`$ as appropriate to the scientific
question.

## Average results across seeds

Rule learning depends on random fold assignments and can vary across
repeated runs. To summarize estimates from several runs of the same
analysis with different random seeds,
[`average_across_seeds()`](https://allicodi.github.io/drotr/reference/average_across_seeds.md)
averages the fold-aggregated estimates and combines the reported
standard errors across runs. In this example, each seed refits the
entire procedure with a new random split, using the same *Shigella*
predictors. This is a descriptive summary of repeated fits, **not** a
substitute for inference over independent datasets or for the
data-adaptive confidence intervals produced by an individual fit.

``` r

Z_shigella <- c("shigella_new", "shigella_bin")

results_list <- vector(mode = "list", length = 3)

for(seed in 1:3){
  
  set.seed(seed)
  
  results_shigella_seed <- estimate_OTR(df = abcd_data,
                                   Y_name = Y, 
                                   A_name = A, 
                                   W_list = W, 
                                   Z_list = Z_shigella, 
                                   id_name = "pid",
                                   sl.library.outcome = sl.library.outcome,
                                   sl.library.treatment = sl.library.treatment,
                                   sl.library.missingness = sl.library.missingness,
                                   sl.library.CATE = sl.library.CATE,
                                   threshold = t,
                                   k_folds = 10,
                                   outcome_type = "gaussian")
  
  results_list[[seed]] <- results_shigella_seed
  
}

average_across_seeds(results_list = results_list, threshold = t)
```

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] drotr_1.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6           knitr_1.52          rlang_1.3.0        
#>  [4] xfun_0.62           otel_0.2.0          textshaping_1.0.5  
#>  [7] SuperLearner_2.0-42 jsonlite_2.0.0      htmltools_0.5.9    
#> [10] ragg_1.5.2          sass_0.4.10         rmarkdown_2.32     
#> [13] evaluate_1.0.5      jquerylib_0.1.4     fastmap_1.2.0      
#> [16] foreach_1.5.2       yaml_2.3.12         lifecycle_1.0.5    
#> [19] compiler_4.6.1      codetools_0.2-20    fs_2.1.0           
#> [22] systemfonts_1.3.2   digest_0.6.39       R6_2.6.1           
#> [25] splines_4.6.1       gam_1.22-7          bslib_0.12.0       
#> [28] tools_4.6.1         nnls_1.6            iterators_1.0.14   
#> [31] pkgdown_2.2.1       cachem_1.1.0        desc_1.4.3
```

Kennedy, Edward H. 2023. *Towards Optimal Doubly Robust Estimation of
Heterogeneous Causal Effects*. arXiv.
<https://doi.org/10.48550/arXiv.2004.14497>.

Luedtke, Alexander R., and Mark J. van der Laan. 2016. “Super-Learning
of an Optimal Dynamic Treatment Rule.” *The International Journal of
Biostatistics* 12 (1): 305. <https://doi.org/10.1515/ijb-2015-0052>.

Montoya, Lina M., Mark J. van der Laan, Alexander R. Luedtke, Jennifer
L. Skeem, Jeremy R. Coyle, and Maya L. Petersen. 2023. “The Optimal
Dynamic Treatment Rule Superlearner: Considerations, Performance, and
Application to Criminal Justice Interventions.” *The International
Journal of Biostatistics* 19 (1): 217–38.
<https://doi.org/10.1515/ijb-2020-0127>.
