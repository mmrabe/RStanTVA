
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
#> Trying to compile a simple C file
#> Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
#> using C compiler: ‘Apple clang version 16.0.0 (clang-1600.0.26.6)’
#> using SDK: ‘MacOSX15.2.sdk’
#> clang -arch arm64 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/Rcpp/include/"  -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/RcppEigen/include/"  -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/RcppEigen/include/unsupported"  -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/BH/include" -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/StanHeaders/include/src/"  -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/StanHeaders/include/"  -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/RcppParallel/include/"  -I"/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c foo.c -o foo.o
#> In file included from <built-in>:1:
#> In file included from /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22:
#> In file included from /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/RcppEigen/include/Eigen/Dense:1:
#> In file included from /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/RcppEigen/include/Eigen/Core:19:
#> /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: 'cmath' file not found
#>   679 | #include <cmath>
#>       |          ^~~~~~~
#> 1 error generated.
#> make: *** [foo.o] Error 1

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
#> 77.5144282  0.5931664 21.9278293 13.1074409
```

Fit `tva_model` to the `tva_data` using maximum-likelihood estimation
(MLE):

``` r
tva_fit <- sampling(tva_model, tva_data)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.007821 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 78.21 seconds.
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
#> Chain 1:  Elapsed Time: 103.694 seconds (Warm-up)
#> Chain 1:                86.576 seconds (Sampling)
#> Chain 1:                190.27 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.006995 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 69.95 seconds.
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
#> Chain 2:  Elapsed Time: 96.243 seconds (Warm-up)
#> Chain 2:                63.855 seconds (Sampling)
#> Chain 2:                160.098 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.007001 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 70.01 seconds.
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
#> Chain 3:  Elapsed Time: 100.979 seconds (Warm-up)
#> Chain 3:                91.461 seconds (Sampling)
#> Chain 3:                192.44 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.006709 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 67.09 seconds.
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
#> Chain 4:  Elapsed Time: 98.614 seconds (Warm-up)
#> Chain 4:                105.161 seconds (Sampling)
#> Chain 4:                203.775 seconds (Total)
#> Chain 4:
tva_fit
#> StanTVA model with 6 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
#> 
#> Model configuration:
#> locations = 4
#> task = "pr"
#> regions = list()
#> C_mode = "equal"
#> w_mode = "locations"
#> t0_mode = "gaussian"
#> K_mode = "free"
#> max_K = 4
#> allow_guessing = FALSE
#> parallel = FALSE
#> save_log_lik = FALSE
#> sanity_checks = FALSE
#> debug_neginf_loglik = FALSE
#> Warning: Unknown or uninitialised column: `param`.
#> 
#> Global parameters:
#>         mean se_mean    sd  2.5%   25%   50%   75%  97.5%   n_eff Rhat
#> C      81.91    0.39 18.61 52.90 68.71 79.25 92.25 124.91 2239.17    1
#> w[1]    0.30    0.00  0.02  0.26  0.28  0.30  0.31   0.34 5021.88    1
#> w[2]    0.29    0.00  0.02  0.25  0.28  0.29  0.31   0.34 4434.08    1
#> w[3]    0.19    0.00  0.02  0.16  0.18  0.19  0.21   0.23 5305.05    1
#> w[4]    0.21    0.00  0.02  0.18  0.20  0.21  0.23   0.25 5375.13    1
#> pK[1]   0.16    0.00  0.02  0.12  0.15  0.16  0.18   0.21 4377.73    1
#> pK[2]   0.18    0.00  0.02  0.14  0.17  0.18  0.20   0.23 4227.23    1
#> pK[3]   0.28    0.00  0.03  0.22  0.26  0.28  0.30   0.34 4741.19    1
#> pK[4]   0.20    0.00  0.03  0.15  0.18  0.20  0.22   0.26 4519.14    1
#> pK[5]   0.17    0.00  0.03  0.12  0.15  0.17  0.19   0.24 3706.31    1
#> mu0    21.65    0.11  4.79 11.59 18.55 21.93 25.08  30.34 1872.49    1
#> sigma0 13.39    0.06  3.14  6.81 11.47 13.47 15.41  19.55 2366.34    1
#> alpha   0.62    0.00  0.15  0.36  0.52  0.61  0.72   0.97 3867.74    1
```
