
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

## Loading the library

``` r
library(RStanTVA)
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

## Simple example

Load data of simulated subject 20 in condition `high` from the parameter
recovery study data set `tva_recovery`:

``` r
tva_data <- tva_recovery %>% filter(subject == 20 & condition == "high")
tva_data
#> # A tibble: 240 × 8
#> # Groups:   subject [1]
#>    subject trial     T condition S[,1]  [,2]  [,3] D[,1] R[,1] true_values$t[,1]
#>      <int> <int> <dbl> <fct>     <int> <int> <int> <int> <int>             <dbl>
#>  1      20     1    50 high          1     1     1     0     0              31.9
#>  2      20     2   200 high          1     1     1     0     0              52.1
#>  3      20     3   200 high          1     1     1     1     0              68.4
#>  4      20     5    50 high          1     1     1     0     0              85.9
#>  5      20     6   200 high          1     1     1     0     1              58.8
#>  6      20     7    50 high          1     1     1     0     0             114. 
#>  7      20    12    10 high          1     1     1     0     0              30.5
#>  8      20    13   200 high          1     1     1     1     0              46.5
#>  9      20    16   200 high          0     1     1     0     0              NA  
#> 10      20    17    50 high          1     1     1     0     1              71.2
#> # ℹ 230 more rows
#> # ℹ 13 more variables: S[4:6] <int>, D[2:6] <int>, R[2:6] <int>,
#> #   true_values$t[2:6] <dbl>, true_values$v <dbl[,6]>, $t0 <dbl>, $K <int>,
#> #   $C <dbl>, $alpha <dbl>, $w <dbl[,4]>, $mu0 <dbl>, $sigma0 <dbl>,
#> #   $pK <dbl[,5]>
```

Generate a StanTVA model for partial report of 6 display locations with
Gaussian $t_0$ and a free $K$ distribution:

``` r
tva_model <- stantva_model(
  locations = 6,
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
#>   - locations = 6
#>   - task = "pr"
#>   - regions = list()
#>   - C_mode = "equal"
#>   - w_mode = "locations"
#>   - t0_mode = "gaussian"
#>   - K_mode = "free"
#>   - max_K = 6
#>   - allow_guessing = FALSE
#>   - parallel = TRUE
#>   - save_log_lik = FALSE
#>   - sanity_checks = FALSE
#>   - debug_neginf_loglik = FALSE
```

The generated Stan code looks like this:

``` stan
/************************************************************************************
 *  StanTVA
 *  =======
 *  This is a StanTVA program, generated with RStanTVA (v0.1.2) Please cite as:
 *  
 *  Rabe M, Kyllingsbæk S (2025). _RStanTVA: TVA Models in 'Stan' using 'R'
 *  and 'StanTVA'_. R package version 0.1.2,
 *  <https://github.com/mmrabe/RStanTVA>.
 *  
 *  Configuration
 *  =============
 *  - formula = NULL
 *  - locations = 6
 *  - task = pr
 *  - regions = list()
 *  - C_mode = equal
 *  - w_mode = locations
 *  - t0_mode = gaussian
 *  - K_mode = free
 *  - max_K = 6
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
        vector[6] s = rep_vector(C, 6);
        array[nS] int Ss = get_matches(S);
        vector[6] w_alpha = w;
        for(i in 1:6) if(D[i]) w_alpha[i] *= alpha;
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
    array[N,6] int<lower=0,upper=1> S;
    array[N,6] int<lower=0,upper=1> R;
    array[N,6] int<lower=0,upper=1> D;
}
transformed data {
    array[N] int nS;
    for(i in 1:N) nS[i] = sum(S[i,]);
    int total_nS = sum(nS);
}
parameters {
    real<lower=machine_precision()> C;
    simplex[6] w;
    simplex[7] pK;
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
#> 45.4031657  0.6393046 12.5131284 12.9087662
```

Fit `tva_model` to the `tva_data` using maximum-likelihood estimation
(MLE):

``` r
tva_fit <- sampling(tva_model, tva_data)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.003363 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 33.63 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 49.438 seconds (Warm-up)
#> Chain 1:                51.695 seconds (Sampling)
#> Chain 1:                101.133 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: Rejecting initial value:
#> Chain 2:   Gradient evaluated at the initial value is not finite.
#> Chain 2:   Stan can't start sampling from this initial value.
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.00333 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 33.3 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 51.803 seconds (Warm-up)
#> Chain 2:                50.808 seconds (Sampling)
#> Chain 2:                102.611 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: Rejecting initial value:
#> Chain 3:   Gradient evaluated at the initial value is not finite.
#> Chain 3:   Stan can't start sampling from this initial value.
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.003326 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 33.26 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 45.793 seconds (Warm-up)
#> Chain 3:                44.434 seconds (Sampling)
#> Chain 3:                90.227 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.005715 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 57.15 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 48.963 seconds (Warm-up)
#> Chain 4:                48.202 seconds (Sampling)
#> Chain 4:                97.165 seconds (Total)
#> Chain 4:
tva_fit
#> StanTVA model with 6 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
#> 
#> Model configuration:
#> locations = 6
#> task = "pr"
#> regions = list()
#> C_mode = "equal"
#> w_mode = "locations"
#> t0_mode = "gaussian"
#> K_mode = "free"
#> max_K = 6
#> allow_guessing = FALSE
#> parallel = TRUE
#> save_log_lik = FALSE
#> sanity_checks = FALSE
#> debug_neginf_loglik = FALSE
#> Warning: Unknown or uninitialised column: `param`.
#> 
#> Global parameters:
#>         mean se_mean    sd  2.5%   25%   50%   75% 97.5%   n_eff Rhat
#> C      49.70    0.23 12.21 32.13 41.06 47.45 55.95 80.02 2875.77    1
#> w[1]    0.18    0.00  0.02  0.14  0.16  0.18  0.20  0.23 5444.68    1
#> w[2]    0.18    0.00  0.02  0.13  0.16  0.18  0.19  0.23 5940.18    1
#> w[3]    0.15    0.00  0.02  0.11  0.13  0.14  0.16  0.19 5947.21    1
#> w[4]    0.19    0.00  0.02  0.14  0.17  0.19  0.20  0.24 5652.26    1
#> w[5]    0.18    0.00  0.02  0.14  0.16  0.18  0.20  0.23 5816.71    1
#> w[6]    0.13    0.00  0.02  0.09  0.11  0.13  0.14  0.17 6775.96    1
#> pK[1]   0.15    0.00  0.03  0.09  0.12  0.14  0.17  0.21 5670.64    1
#> pK[2]   0.10    0.00  0.02  0.06  0.08  0.09  0.11  0.15 6611.23    1
#> pK[3]   0.18    0.00  0.04  0.11  0.15  0.18  0.20  0.26 5295.95    1
#> pK[4]   0.20    0.00  0.04  0.12  0.17  0.20  0.22  0.28 5309.87    1
#> pK[5]   0.14    0.00  0.04  0.08  0.11  0.14  0.16  0.22 5307.71    1
#> pK[6]   0.12    0.00  0.03  0.07  0.10  0.12  0.14  0.20 4979.90    1
#> pK[7]   0.12    0.00  0.03  0.07  0.10  0.12  0.14  0.19 5219.54    1
#> mu0    13.99    0.13  6.82  1.92  8.81 13.83 18.83 27.40 2612.09    1
#> sigma0 15.24    0.12  6.16  3.56 11.01 15.25 19.27 27.41 2742.11    1
#> alpha   0.69    0.00  0.18  0.39  0.56  0.67  0.80  1.08 4607.48    1
```

## Nested example

Load data of simulated subject 20 in both condition (`high` and `low`)
from the parameter recovery study data set `tva_recovery`:

``` r
tva_data_nested <- tva_recovery %>% filter(subject == 20)
tva_data_nested
#> # A tibble: 480 × 8
#> # Groups:   subject [1]
#>    subject trial     T condition S[,1]  [,2]  [,3] D[,1] R[,1] true_values$t[,1]
#>      <int> <int> <dbl> <fct>     <int> <int> <int> <int> <int>             <dbl>
#>  1      20     1    50 high          1     1     1     0     0              31.9
#>  2      20     2   200 high          1     1     1     0     0              52.1
#>  3      20     3   200 high          1     1     1     1     0              68.4
#>  4      20     4    50 low           1     1     1     1     0              34.5
#>  5      20     5    50 high          1     1     1     0     0              85.9
#>  6      20     6   200 high          1     1     1     0     1              58.8
#>  7      20     7    50 high          1     1     1     0     0             114. 
#>  8      20     8   200 low           1     1     1     0     1              81.9
#>  9      20     9    50 low           1     1     1     1     0             207. 
#> 10      20    10   100 low           1     1     1     0     1              39.3
#> # ℹ 470 more rows
#> # ℹ 13 more variables: S[4:6] <int>, D[2:6] <int>, R[2:6] <int>,
#> #   true_values$t[2:6] <dbl>, true_values$v <dbl[,6]>, $t0 <dbl>, $K <int>,
#> #   $C <dbl>, $alpha <dbl>, $w <dbl[,4]>, $mu0 <dbl>, $sigma0 <dbl>,
#> #   $pK <dbl[,5]>
```

Generate a StanTVA model for partial report of 6 display locations with
Gaussian $t_0$ and a free $K$ distribution:

``` r
tva_model_nested <- stantva_model(
  formula = list(
    log(C) ~ 1 + condition
  ),
  locations = 6,
  task = "pr",
  w_mode = "locations",
  t0_mode = "gaussian",
  K_mode = "free",
  sanity_checks = FALSE
)
```

Note that we are allowing $C$ to vary between experimental conditions.
This will add another “layer” to the model, which implements $C$ in
trial $i$, $C_i$, as:

$$
\log\left(C_i\right) = \alpha_C+\beta_CX_i
$$

As a consequence, $C$ is no longer estimated as a single invariant
parameter but as the exp-sum of an intercept $\alpha_C$, and the product
of slope $\beta_C$ and experimental condition $X_i$, coded here using
standard treatment contrasts.

The generated Stan code looks like this:

``` stan
/************************************************************************************
 *  StanTVA
 *  =======
 *  This is a StanTVA program, generated with RStanTVA (v0.1.2) Please cite as:
 *  
 *  Rabe M, Kyllingsbæk S (2025). _RStanTVA: TVA Models in 'Stan' using 'R'
 *  and 'StanTVA'_. R package version 0.1.2,
 *  <https://github.com/mmrabe/RStanTVA>.
 *  
 *  Configuration
 *  =============
 *  - formula = list(log(C) ~ 1 + condition)
 *  - locations = 6
 *  - task = pr
 *  - regions = list()
 *  - C_mode = equal
 *  - w_mode = locations
 *  - t0_mode = gaussian
 *  - K_mode = free
 *  - max_K = 6
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
    vector calculate_v(data int nS, data array[] int S, data array[] int D, vector w, real alpha, real C) {
        vector[6] s = rep_vector(C, 6);
        array[nS] int Ss = get_matches(S);
        vector[6] w_alpha = w;
        for(i in 1:6) if(D[i]) w_alpha[i] *= alpha;
        vector[nS] v = s[Ss] .* w_alpha[Ss] / sum(w_alpha[Ss]);
        for(i in 1:nS) if(v[i] < machine_precision()) v[i] = machine_precision();
        return v/1000.0;
    }
    real log_lik_single(data array[] int S, data array[] int D, data array[] int R, data int nS, data real T, vector w, real alpha, vector pK, real mu0, real sigma0, real C) {
        real log_lik;
        vector[nS] v = calculate_v(nS, S, D, to_vector(w), alpha, C);
        log_lik = tva_pr_log(R, S, D, T, [mu0, sigma0]', pK, v);
        return log_lik;
    }
}
data {
    int<lower=1> N;
    array[N] real<lower=0> T;
    array[N,6] int<lower=0,upper=1> S;
    array[N,6] int<lower=0,upper=1> R;
    array[N,6] int<lower=0,upper=1> D;
    int<lower=0> M_C;
    int<lower=0,upper=M_C> int_C;
    matrix[N,M_C] X;
    array[M_C] int map_C;
}
transformed data {
    array[N] int nS;
    for(i in 1:N) nS[i] = sum(S[i,]);
    int total_nS = sum(nS);
    int M = M_C;
}
parameters {
    simplex[6] w;
    simplex[7] pK;
    real mu0;
    real<lower=machine_precision()> sigma0;
    real<lower=machine_precision()> alpha;
    vector[M] b;
}
transformed parameters {
    vector<lower=machine_precision()>[N] C;
    {
        C = X[,map_C] * b[map_C];
        C = exp(C);
    }
    for(i in 1:N) if(is_nan(C[i]) || is_inf(C[i])) reject("Rejecting proposal because C[",i,"] = ",C[i]," !");
}
model {
    if(int_C) {
        b[map_C[:(int_C-1)]] ~ normal(0.0,5.0);
        b[map_C[(int_C+1):]] ~ normal(0.0,5.0);
        b[map_C[int_C]] ~ normal(0.0,10.0);
    } else {
        b[map_C] ~ normal(0.0,5.0);
    }
    w ~ lognormal(0, 0.5);
    pK ~ lognormal(0, 0.5);
    mu0 ~ normal(20, 15);
    sigma0 ~ gamma(2, 0.2);
    alpha ~ lognormal(-0.4, 0.6);
    // likelihood (only if prior != 0)
    if(target() != negative_infinity()) {
        for(i in 1:N) target += log_lik_single(S[i], D[i], R[i], nS[i], T[i], to_vector(w), alpha, to_vector(pK), mu0, sigma0, C[i]);
    }
}
```

Fit `tva_model_nested` to the `tva_data_nested` using maximum-likelihood
estimation (MLE):

``` r
tva_fit_nested <- sampling(tva_model_nested, tva_data_nested)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.006813 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 68.13 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 96.61 seconds (Warm-up)
#> Chain 1:                95.181 seconds (Sampling)
#> Chain 1:                191.791 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: Rejecting initial value:
#> Chain 2:   Gradient evaluated at the initial value is not finite.
#> Chain 2:   Stan can't start sampling from this initial value.
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.006249 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 62.49 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 99.622 seconds (Warm-up)
#> Chain 2:                96.655 seconds (Sampling)
#> Chain 2:                196.277 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: Rejecting initial value:
#> Chain 3:   Gradient evaluated at the initial value is not finite.
#> Chain 3:   Stan can't start sampling from this initial value.
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.005611 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 56.11 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 105.754 seconds (Warm-up)
#> Chain 3:                91.585 seconds (Sampling)
#> Chain 3:                197.339 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.006701 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 67.01 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 107.096 seconds (Warm-up)
#> Chain 4:                89.37 seconds (Sampling)
#> Chain 4:                196.466 seconds (Total)
#> Chain 4:
tva_fit_nested
#> StanTVA model with 6 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
#> 
#> Model configuration:
#> formula = list(log(C) ~ 1 + condition)
#> locations = 6
#> task = "pr"
#> regions = list()
#> C_mode = "equal"
#> w_mode = "locations"
#> t0_mode = "gaussian"
#> K_mode = "free"
#> max_K = 6
#> allow_guessing = FALSE
#> parallel = TRUE
#> save_log_lik = FALSE
#> sanity_checks = FALSE
#> debug_neginf_loglik = FALSE
#> 
#> Global parameters:
#>         mean se_mean   sd 2.5%   25%   50%   75% 97.5%   n_eff Rhat
#> w[1]    0.18    0.00 0.02 0.15  0.17  0.18  0.19  0.22 6518.94    1
#> w[2]    0.17    0.00 0.02 0.13  0.15  0.17  0.18  0.20 6984.50    1
#> w[3]    0.15    0.00 0.02 0.12  0.13  0.14  0.16  0.18 7359.58    1
#> w[4]    0.17    0.00 0.02 0.14  0.16  0.17  0.19  0.21 7534.26    1
#> w[5]    0.20    0.00 0.02 0.17  0.19  0.20  0.21  0.24 6847.21    1
#> w[6]    0.13    0.00 0.01 0.10  0.12  0.13  0.14  0.16 6935.49    1
#> pK[1]   0.13    0.00 0.02 0.09  0.12  0.13  0.15  0.18 6200.28    1
#> pK[2]   0.10    0.00 0.02 0.06  0.08  0.10  0.11  0.14 6772.40    1
#> pK[3]   0.20    0.00 0.04 0.13  0.18  0.20  0.23  0.28 4815.87    1
#> pK[4]   0.20    0.00 0.04 0.13  0.18  0.20  0.23  0.29 7055.99    1
#> pK[5]   0.14    0.00 0.03 0.08  0.11  0.14  0.16  0.22 6319.33    1
#> pK[6]   0.11    0.00 0.03 0.06  0.09  0.11  0.13  0.18 5601.24    1
#> pK[7]   0.11    0.00 0.03 0.06  0.09  0.11  0.13  0.18 5376.86    1
#> mu0    14.19    0.10 5.29 5.20 10.06 13.98 17.72 24.97 2919.25    1
#> sigma0 12.12    0.08 4.66 3.14  8.88 12.16 15.55 21.00 3464.74    1
#> alpha   0.74    0.00 0.14 0.50  0.64  0.73  0.82  1.02 7651.98    1
#> Warning: Unknown or uninitialised column: `param`.
#> 
#> Fixed effects:
#>                 mean se_mean   sd 2.5%  25%  50%  75% 97.5%   n_eff Rhat
#> C_Intercept     3.50       0 0.17 3.20 3.38 3.49 3.61  3.85 3192.88    1
#> C_conditionhigh 0.39       0 0.17 0.06 0.28 0.38 0.50  0.73 6648.00    1
```

## Hierarchical example

Generate a hierarchical TVA model:

``` r


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
    log(C) ~ 1 + condition + (1 + condition | C_alpha | subject),
    log(w) ~ 1 + (1 | subject),
    log(pK) ~ 1 + (1 | subject),
    mu0 ~ 1 + (1 | subject),
    log(sigma0) ~ 1 + (1 | subject),
    log(alpha) ~ 1 + (1 | C_alpha | subject)
  ),
  locations = 6,
  task = "pr",
  w_mode = "locations",
  t0_mode = "gaussian",
  K_mode = "free",
  priors = priors
)
```

Fit hierarchical `tva_hierarchical_model` to the first 200 trials of the
first 10 subjects of the `tva_recovery` data set:

``` r

tva_hierarchical_subset <- tva_recovery %>% filter(subject <= 10 & trial <= 200)

tva_hierarchical_fit <- sampling(tva_hierarchical_model, tva_hierarchical_subset, refresh = 20)

tva_hierarchical_fit
```
