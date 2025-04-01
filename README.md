
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RStanTVA

<!-- badges: start -->

[![R](https://github.com/mmrabe/RStanTVA/actions/workflows/rpkg.yml/badge.svg)](https://github.com/mmrabe/RStanTVA/actions/workflows/rpkg.yml)
<!-- badges: end -->

RStanTVA is an R package containing the StanTVA library and numerous
convenience functions for generating, compiling, fitting, and analyzing
(Stan)TVA models.

## Installation

You can install the development version of RStanTVA from
[GitHub](https://github.com/mmrabe/RStanTVA) with:

``` r
remotes::install_github("mmrabe/RStanTVA")
```

## Example

Load the R package:

``` r
library(RStanTVA)
#> Loading required package: rstan
#> Loading required package: StanHeaders
#> 
#> rstan version 2.32.6 (Stan version 2.32.2)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
#> For within-chain threading using `reduce_sum()` or `map_rect()` Stan functions,
#> change `threads_per_chain` option:
#> rstan_options(threads_per_chain = 1)
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ tidyr::extract() masks rstan::extract()
#> ✖ dplyr::filter()  masks stats::filter()
#> ✖ dplyr::lag()     masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

Load example data from the parameter recovery study:

``` r
tva_data <- tva_recovery %>% filter(subject == 20)
tva_data
#> # A tibble: 512 × 9
#> # Groups:   subject [1]
#>    subject trial     T condition type  S[,1]  [,2] D[,1] R[,1] true_values$t[,1]
#>      <int> <int> <dbl> <fct>     <fct> <int> <int> <int> <int>             <dbl>
#>  1      20     1    10 high      WR        1     1     0     0              26.1
#>  2      20     2    10 low       WR        1     1     0     0              56.4
#>  3      20     3    50 high      WR        1     1     0     0             139. 
#>  4      20     4    50 high      PR        1     1     0     0              72.6
#>  5      20     5    50 low       WR        1     1     0     0              62.9
#>  6      20     6    50 low       PR        1     1     0     0              95.7
#>  7      20     7   100 high      WR        1     1     0     1              37.8
#>  8      20     8   100 high      PR        1     1     1     0             118. 
#>  9      20     9   100 low       WR        1     1     0     0             125. 
#> 10      20    10   100 low       PR        1     1     1     0              77.3
#> # ℹ 502 more rows
#> # ℹ 13 more variables: S[3:4] <int>, D[2:4] <int>, R[2:4] <int>,
#> #   true_values$t[2:6] <dbl>, true_values$v <dbl[,6]>, $t0 <dbl>, $K <int>,
#> #   $C <dbl>, $alpha <dbl>, $w <dbl[,4]>, $mu0 <dbl>, $sigma0 <dbl>,
#> #   $pK <dbl[,5]>
```

Generate a StanTVA model for partial report of 6 display locations with
Gaussian $t_0$ and a free $K$ distribution:

``` r
library(RStanTVA)

tva_model <- stantva_model(
  location = 4,
  task = "pr",
  w_mode = "locations",
  t0_mode = "gaussian",
  K_mode = "free",
  sanity_checks = FALSE
)
```

``` r
tva_model
#> StanTVA model with 6 free parameter(s) and the following configuration:
#>   - locations = 4
#>   - task = "pr"
#>   - regions = list()
#>   - C_mode = "equal"
#>   - w_mode = "locations"
#>   - t0_mode = "gaussian"
#>   - K_mode = "free"
#>   - max_K = 4
#>   - allow_guessing = FALSE
#>   - parallel = FALSE
#>   - save_log_lik = FALSE
#>   - sanity_checks = FALSE
#>   - debug_neginf_loglik = FALSE
```

The generated Stan model looks like this:

``` stan
/************************************************************************************
 *  StanTVA
 *  =======
 *  This is a StanTVA program, generated with RStanTVA (v0.1.2) Please cite as:
 *  
 *  Rabe M, Kyllingsbæk S (2025). _RStanTVA: TVA Models in Stan using R and
 *  StanTVA_. R package version 0.1.2.
 *  
 *  Configuration
 *  =============
 *  - formula = NULL
 *  - locations = 4
 *  - task = pr
 *  - regions = list()
 *  - C_mode = equal
 *  - w_mode = locations
 *  - t0_mode = gaussian
 *  - K_mode = free
 *  - max_K = 4
 *  - allow_guessing = FALSE
 *  - parallel = FALSE
 *  - save_log_lik = FALSE
 *  - priors = NULL
 *  - sanity_checks = FALSE
 *  - debug_neginf_loglik = FALSE
 *  
 *  License
 *  =======
 *  This program is licensed under the GNU General Public License 3. For a copy of
 *  the license agreement, see: https://www.gnu.org/licenses/gpl-3.0.html
 ************************************************************************************/

functions {
    #include tva.stan
    #include freeK.stan
    #include gaussiant0.stan
    vector calculate_v(data int nS, data array[] int S, data array[] int D, vector w, real C, real alpha) {
        vector[4] s = rep_vector(C, 4);
        array[nS] int Ss = get_matches(S);
        vector[4] w_alpha = w;
        for(i in 1:4) if(D[i]) w_alpha[i] *= alpha;
        vector[nS] v = s[Ss] .* w_alpha[Ss] / sum(w_alpha[Ss]);
        for(i in 1:nS) if(v[i] < machine_precision()) v[i] = machine_precision();
        return v/1000.0;
    }
    real log_lik_single(data array[] int S, data array[] int D, data array[] int R, data int nS, data real T, vector w, real C, real alpha, vector pK, real mu0, real sigma0) {
        real log_lik;
        vector[nS] v = calculate_v(nS, S, D, to_vector(w), C, alpha);
        log_lik = tva_pr_log(R, S, D, T, [mu0, sigma0]', pK, v);
        return log_lik;
    }
}
data {
    int<lower=1> N;
    array[N] real<lower=0> T;
    array[N,4] int<lower=0,upper=1> S;
    array[N,4] int<lower=0,upper=1> R;
    array[N,4] int<lower=0,upper=1> D;
}
transformed data {
    array[N] int nS;
    for(i in 1:N) nS[i] = sum(S[i,]);
    int total_nS = sum(nS);
}
parameters {
    real<lower=machine_precision()> C;
    simplex[4] w;
    simplex[5] pK;
    real mu0;
    real<lower=machine_precision()> sigma0;
    real<lower=machine_precision()> alpha;
}
model {
    C ~ gamma(3.5, 0.035);
    w ~ lognormal(0, 0.5);
    pK ~ lognormal(0, 0.5);
    mu0 ~ normal(20, 15);
    sigma0 ~ gamma(2, 0.2);
    alpha ~ lognormal(-0.4, 0.6);
    // likelihood (only if prior != 0)
    if(target() != negative_infinity()) {
        for(i in 1:N) target += log_lik_single(S[i], D[i], R[i], nS[i], T[i], to_vector(w), C, alpha, to_vector(pK), mu0, sigma0);
    }
}
```

Fit `tva_model` to the `tva_data` using maximum-likelihood estimation
(MLE):

``` r
tva_fit_mle <- optimizing(tva_model, tva_data)
tva_fit_mle$par[c("C","alpha","mu0","sigma0")]
#>          C      alpha        mu0     sigma0 
#> 77.5300589  0.5932424 21.9315554 13.1077232
```

Fit `tva_model` to the `tva_data` using maximum-likelihood estimation
(MLE):

``` r
tva_fit <- sampling(tva_model, tva_data)
tva_fit
```

Generate a hierarchical TVA model:

``` r

rstan_options(threads_per_chain = parallel::detectCores())


priors <-
  prior(normal(0,.07),dpar=C)+
  prior(normal(4,.2),dpar=C,coef=Intercept)+
  prior(normal(0,.07),dpar=alpha)+
  prior(normal(-0.2,.1),dpar=alpha,coef=Intercept)+
  prior(normal(0,.03),dpar=pK)+
  prior(normal(0,.1),dpar=pK,coef=Intercept)+
  prior(normal(0,5),dpar=mu0)+
  prior(normal(30,15),dpar=mu0,coef=Intercept)+
  prior(normal(0,.04),dpar=sigma0)+
  prior(normal(0,.2),dpar=sigma0,coef=Intercept)+
  prior(normal(0,.05),dpar=w)+
  prior(normal(0,0.1),dpar=w,coef=Intercept)

tva_hierarchical_model <- stantva_model(
  formula = list(
    log(C) ~ 1 + (1 | C_alpha | subject),
    log(w) ~ 1 + (1 | subject),
    log(pK) ~ 1 + (1 | subject),
    mu0 ~ 1 + (1 | subject),
    log(sigma0) ~ 1 + (1 | subject),
    log(alpha) ~ 1 + (1 | C_alpha | subject)
  ),
  location = 4,
  task = "pr",
  w_mode = "locations",
  t0_mode = "gaussian",
  K_mode = "free",
  priors = priors,
  sanity_checks = FALSE,
  parallel = TRUE
)
```

Fit hierarchical `tva_hierarchical_model` to the first 200 trials of the
first 10 subjects of `tva_recovery`:

``` r

tva_hierarchical_subset <- tva_recovery %>% filter(subject <= 10 & trial <= 200)

tva_hierarchical_fit <- sampling(tva_hierarchical_model, tva_hierarchical_subset, cores = 4, chains = 4, refresh = 20)

tva_hierarchical_fit
```

    #> StanTVA model with 6 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
    #> 
    #> Model configuration:
    #> formula = list(log(C) ~ 1 + (1 | C_alpha | subject), log(w) ~ 1 + (1 | subject), log(pK) ~ 1 + (1 | subject), mu0 ~ 1 + (1 | subject), log(sigma0) ~ 1 + (1 | subject), log(alpha) ~ 1 + (1 | C_alpha | subject))
    #> locations = 4
    #> task = "pr"
    #> regions = list()
    #> C_mode = "equal"
    #> w_mode = "locations"
    #> t0_mode = "gaussian"
    #> K_mode = "free"
    #> max_K = 4
    #> allow_guessing = FALSE
    #> parallel = TRUE
    #> save_log_lik = FALSE
    #> priors = structure(list(prior = c("normal(0, 0.07)", "normal(4, 0.2)", "normal(0, 0.07)", "normal(-0.2, 0.1)", "normal(0, 0.03)", "normal(0, 0.1)", "normal(0, 5)", "normal(30, 15)", "normal(0, 0.04)", "normal(0, 0.2)", "normal(0, 0.05)", "normal(0, 0.1)"), class = c("b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b"), coef = c("", "Intercept", "", "Intercept", "", "Intercept", "", "Intercept", "", "Intercept", "", "Intercept"), group = c("", "", "", "", "", "", "", "", "", "", "", ""), resp = c("",  "", "", "", "", "", "", "", "", "", "", ""), dpar = c("C", "C", "alpha", "alpha", "pK", "pK", "mu0", "mu0", "sigma0", "sigma0", "w", "w"), nlpar = c("", "", "", "", "", "", "", "", "", "", "", ""), lb = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_), ub = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,  NA_character_, NA_character_, NA_character_, NA_character_), source = c("user", "user", "user", "user", "user", "user", "user", "user", "user", "user", "user", "user")), row.names = c(NA, -12L), class = c("brmsprior", "data.frame"))
    #> sanity_checks = FALSE
    #> debug_neginf_loglik = FALSE
    #> 
    #> Fixed effects:
    #>                   mean se_mean   sd  2.5%   25%   50%   75% 97.5%   n_eff Rhat
    #> C_Intercept       3.97    0.00 0.08  3.82  3.92  3.96  4.01  4.12 2108.39    1
    #> w_Intercept[1]    0.14    0.00 0.05  0.03  0.10  0.14  0.18  0.25 2437.63    1
    #> w_Intercept[2]    0.07    0.00 0.06 -0.04  0.03  0.07  0.11  0.18 4018.12    1
    #> w_Intercept[3]    0.01    0.00 0.06 -0.11 -0.03  0.00  0.05  0.13 1987.49    1
    #> pK_Intercept[1]  -0.10    0.00 0.07 -0.25 -0.15 -0.10 -0.05  0.04 3207.99    1
    #> pK_Intercept[2]  -0.22    0.00 0.08 -0.37 -0.27 -0.22 -0.17 -0.06 2705.03    1
    #> pK_Intercept[3]   0.32    0.00 0.09  0.14  0.26  0.32  0.39  0.50 1246.13    1
    #> pK_Intercept[4]   0.06    0.00 0.08 -0.10  0.01  0.07  0.12  0.23 2855.05    1
    #> mu0_Intercept     8.50    0.01 0.53  7.43  8.16  8.50  8.84  9.54 2262.46    1
    #> sigma0_Intercept  0.15    0.00 0.22 -0.26  0.00  0.16  0.30  0.57 3568.73    1
    #> alpha_Intercept  -0.35    0.00 0.08 -0.50 -0.39 -0.35 -0.30 -0.20 4555.10    1
    #> 
    #> Hyperparameters on random effects (subject level, N = 10):
    #>                                                  mean se_mean   sd  2.5%   25%
    #> sd(C_subject_Intercept)                          0.17    0.00 0.05  0.06  0.14
    #> sd(alpha_subject_Intercept)                      0.08    0.00 0.06  0.01  0.03
    #> sd(w_subject_Intercept[1])                       0.06    0.00 0.05  0.01  0.03
    #> sd(w_subject_Intercept[2])                       0.09    0.00 0.06  0.01  0.04
    #> sd(w_subject_Intercept[3])                       0.14    0.00 0.06  0.02  0.09
    #> sd(pK_subject_Intercept[1])                      0.11    0.00 0.07  0.01  0.06
    #> sd(pK_subject_Intercept[2])                      0.09    0.00 0.06  0.01  0.05
    #> sd(pK_subject_Intercept[3])                      0.21    0.00 0.09  0.04  0.15
    #> sd(pK_subject_Intercept[4])                      0.08    0.00 0.06  0.01  0.04
    #> sd(mu0_subject_Intercept)                        0.09    0.01 0.07  0.01  0.04
    #> sd(sigma0_subject_Intercept)                     0.11    0.01 0.07  0.01  0.05
    #> cor(C_subject_Intercept,alpha_subject_Intercept) 0.01    0.01 0.37 -0.68 -0.25
    #>                                                   50%  75% 97.5%   n_eff Rhat
    #> sd(C_subject_Intercept)                          0.17 0.20  0.27 1263.38 1.00
    #> sd(alpha_subject_Intercept)                      0.07 0.11  0.22  175.05 1.03
    #> sd(w_subject_Intercept[1])                       0.05 0.09  0.18  174.37 1.03
    #> sd(w_subject_Intercept[2])                       0.08 0.13  0.21  247.85 1.02
    #> sd(w_subject_Intercept[3])                       0.13 0.18  0.26  422.90 1.01
    #> sd(pK_subject_Intercept[1])                      0.10 0.16  0.26  286.77 1.01
    #> sd(pK_subject_Intercept[2])                      0.08 0.13  0.23  234.36 1.02
    #> sd(pK_subject_Intercept[3])                      0.22 0.27  0.38  433.53 1.01
    #> sd(pK_subject_Intercept[4])                      0.07 0.12  0.21  206.43 1.02
    #> sd(mu0_subject_Intercept)                        0.08 0.13  0.25  120.92 1.06
    #> sd(sigma0_subject_Intercept)                     0.10 0.15  0.29  169.21 1.04
    #> cor(C_subject_Intercept,alpha_subject_Intercept) 0.02 0.28  0.70 1906.34 1.00
    #> Warning in .local(x, ...): Model did not converge (Rhat >= 1.05) for 1
    #> parameter(s): sd(mu0_subject_Intercept)
