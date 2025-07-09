
# BAR12: Bayesian Autoregressive Phase 1-2 Design for Cell Therapy  Trials with  Manufacturing Changes

<!-- badges: start -->
<!-- badges: end -->

The goal of **BAR12 design** is to enable seamless dose-finding in early-phase cell therapy trials that incorporate mid-trial manufacturing modifications—or *tweaks*—such as changes to cell culture duration, cytokine cocktails, or donor criteria. 
BAR12 addresses this by using a Bayesian autoregressive model with spike-and-slab priors to distinguish between pre- and post-tweak outcome distributions, allowing dynamic borrowing of information while accounting for possible shifts.

## 📦 Installation

You can install the development version of BAR12:

``` r
remotes::install_github("cyang728/BAR12")
```

## 🔧 Functions

### `BAR12_design()`

### Inputs

- **dosages:** Numeric vector. The standardized dosages being investigated in the trial.
- **p_T_sim:** List of numeric vectors. The true toxicity probabilities for each dose at each stage, used for outcome simulation.
- **p_E_sim:** List of numeric vectors. The true efficacy probabilities for each dose at each stage, used for outcome simulation.
- **num_stages:** Integer. The total number of stages in the trial design.
- **ncohort:** Numeric vector. The number of cohorts to be enrolled in each stage.
- **cohortsize:** Numeric vector. The number of patients per cohort in each stage.
- **startdose:** Integer. The index of the dose level for the first cohort.
- **target_tox:** Numeric. The target toxicity rate, used for decision-making.
- **target_eff:** Numeric. The target efficacy rate.
- **cutoff_tox:** Numeric. The safety cutoff; the posterior probability of toxicity exceeding `target_tox` must not be greater than this value.
- **cutoff_eff:** Numeric. The futility cutoff; the posterior probability of efficacy being less than `target_eff` must not be greater than this value.
- **w00:** Numeric. Utility weight for the outcome of no toxicity and no efficacy.
- **w01:** Numeric. Utility weight for the outcome of no toxicity and efficacy.
- **w10:** Numeric. Utility weight for the outcome of toxicity and no efficacy.
- **w11:** Numeric. Utility weight for the outcome of toxicity and efficacy.
- **max_allocate_dose:** Integer. The maximum number of patients that can be assigned to a single dose level.
- **n_mc_epsilon:** Integer. The number of Monte Carlo samples used to integrate over the random effect parameter.
- **kappa:** Numeric. A randomization parameter controlling the trade-off between exploration and exploitation. Higher values make allocation favor doses with higher utility.
- **L:** Integer. The number of best-performing doses to which patients can be randomized in later stages.
- **seed:** Integer. A random seed for ensuring the reproducibility of simulations.
- **print.out:** Logical. If `TRUE`, prints live updates for each cohort.

### Outputs

A list with the following elements:
- **data:** A data frame containing the complete trial history and final posterior estimates for each dose level.
- **best_dose:** The index of the final recommended optimal dose based on the highest utility among all admissible doses. Returns `0` if no dose is selected.


## 🚀 Example

The following example demonstrates how to use `BAR12_design()` to simulate a BAR12 trial.

``` r
library(BAR12)

out <- BAR12_design(dosages = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                    p_T_sim = list(c(0.01,0.02,0.03,0.03,0.06,0.08),
                                   c(0.01,0.02,0.03,0.05,0.08,0.28)),
                    p_E_sim = list(c(0.02,0.03,0.12,0.19,0.41,0.56),
                                   c(0.06,0.11,0.15,0.21,0.58,0.58)),
                    num_stages = 2,
                    ncohort = c(6, 8), cohortsize = c(3, 3), startdose = 1,
                    target_tox = 0.3, target_eff = 0.2,
                    cutoff_tox = 0.60, cutoff_eff = 0.85,
                    w00 = 40, w01 = 100, w10 = 0, w11 = 60,
                    max_allocate_dose = 100,
                    n_mc_epsilon = 1000,
                    kappa = 1.0,                     
                    L = 2,
                    seed = 1,
                    print.out = FALSE)

out$best_dose   # index of the optimal dose
head(out$data)  # summary data
```

