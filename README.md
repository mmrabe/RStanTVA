
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

## Example data set

``` r
data("tva_recovery")
data("tva_recovery_true_params")
```

`RStanTVA` comes with a CombiTVA data set `tva_recovery` of 50 simulated
subjects with 234 trials each. For each subject and trial, the true
underlying parameters are known exactly, which makes it very useful for
demonstrating the functionality of this package.

Of the 234 trials, half were carried out under the `high` condition, and
the other half under the `low` condition. Imagine that those are maybe
different levels of luminance, some sort of cognitive load or so. The
concrete manipulation is not really relevant but we *do* know that on
average, $C_\mathrm{high}=100.0>C_\mathrm{low}=80.0$, and that across
subjects, $\mathrm{cor}\left(C,\alpha\right)=-0.3$.

## Simple model

Like predecessor software, `RStanTVA` can be used to fit single
parameters to data sets. For example, if we wanted to estimate TVA
parameters to all trials of subject 20 in condition `high`, then we
could first retrieve that subset from the full data set `tva_recovery`:

``` r
tva_data <- tva_recovery %>% filter(subject == 20 & condition == "high")
tva_data
#> # A tibble: 117 × 9
#> # Groups:   subject [1]
#>    subject trial     T condition    nD S[,1]  [,2] D[,1] R[,1] true_values$t[,1]
#>      <int> <int> <dbl> <fct>     <int> <int> <int> <int> <int>             <dbl>
#>  1      20     2  33.3 high          0     1     1     0     1              34.5
#>  2      20     3  16.6 high          0     1     1     0     0              41.4
#>  3      20     4 150   high          4     1     1     1     0              56.8
#>  4      20    11 150   high          3     1     1     0     1              33.2
#>  5      20    12  33.3 high          0     1     1     0     0              22.5
#>  6      20    14  33.3 high          0     1     1     0     0             122. 
#>  7      20    15  50   high          0     1     1     0     0              79.8
#>  8      20    17  33.3 high          0     1     1     0     0              28.1
#>  9      20    18 150   high          3     1     1     0     1              49.1
#> 10      20    19 150   high          4     1     1     0     0              66.1
#> # ℹ 107 more rows
#> # ℹ 13 more variables: S[3:6] <int>, D[2:6] <int>, R[2:6] <int>,
#> #   true_values$t[2:6] <dbl>, true_values$v <dbl[,6]>, $t0 <dbl>, $K <int>,
#> #   $C <dbl>, $alpha <dbl>, $w <dbl[,6]>, $mu0 <dbl>, $sigma0 <dbl>,
#> #   $pK <dbl[,7]>
```

Then we can define the TVA model that we want to fit to the data. Here,
we are dealing with a CombiTVA paradigm, which is effectively a set of
different partial-report configurations. Let us assume a TVA model with
6 display locations, Gaussian $t_0$, and a free $K$ distribution, then
we can generate it as follows:

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

If we wanted to look at the generated Stan code, this would be it:

``` stan
/************************************************************************************
 *  StanTVA
 *  =======
 *  This is a StanTVA program, generated with RStanTVA (v0.2.0) Please cite as:
 *  
 *  Rabe M, Kyllingsbæk S (2025). _RStanTVA: TVA Models in 'Stan' using 'R'
 *  and 'StanTVA'_. R package version 0.2.0,
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

Note that the file will include `tva.stan`, `freeK.stan`, and
`gaussiant0.stan`, all of which are contained in the StanTVA library
embedded in the package. While `tva.stan` contains a number of more
general functions for the likelihood computation for whole and partial
report paradigms, `freeK.stan` and `gaussiant0.stan` are modules that
specify the free $K$ distribution and the Gaussian $t_0$ distribution,
respectively. There will always be exactly three files included. While
`tva.stan` is always required, the other two differ depending on the
selected assumptions for $K$ and $t_0$, respectively. The includes are
located at `stantva_path()` and the RStanTVA function `stantva_model()`
will automatically include those when compiling.

We can now fit `tva_model` to the `tva_data` using maximum-likelihood
estimation (MLE):

``` r
tva_fit_mle <- optimizing(tva_model, tva_data)
tva_fit_mle$par[c("C","alpha","mu0","sigma0")]
#>           C       alpha         mu0      sigma0 
#> 103.6942724   0.9400363  13.9074017   4.4078324
```

… or using Bayesian HMC sampling:

``` r
tva_fit <- sampling(tva_model, tva_data)
```

``` r
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
#>          mean se_mean    sd  2.5%   25%    50%    75%  97.5%   n_eff Rhat
#> C      106.86    0.56 29.58 62.48 85.82 102.58 122.46 176.85 2759.20    1
#> w[1]     0.14    0.00  0.02  0.10  0.13   0.14   0.15   0.18 4640.87    1
#> w[2]     0.22    0.00  0.03  0.17  0.20   0.22   0.24   0.28 4516.45    1
#> w[3]     0.12    0.00  0.02  0.08  0.11   0.12   0.13   0.16 5781.36    1
#> w[4]     0.16    0.00  0.02  0.12  0.14   0.16   0.18   0.21 4003.64    1
#> w[5]     0.19    0.00  0.03  0.14  0.17   0.19   0.21   0.24 4137.62    1
#> w[6]     0.17    0.00  0.02  0.12  0.15   0.17   0.18   0.22 4426.74    1
#> pK[1]    0.11    0.00  0.02  0.07  0.09   0.10   0.12   0.15 4288.52    1
#> pK[2]    0.11    0.00  0.03  0.06  0.09   0.11   0.13   0.17 4105.34    1
#> pK[3]    0.20    0.00  0.04  0.13  0.18   0.20   0.23   0.28 4097.63    1
#> pK[4]    0.22    0.00  0.04  0.15  0.19   0.22   0.24   0.30 4694.32    1
#> pK[5]    0.12    0.00  0.03  0.07  0.10   0.12   0.14   0.18 4192.42    1
#> pK[6]    0.13    0.00  0.03  0.08  0.10   0.12   0.15   0.19 4244.12    1
#> pK[7]    0.12    0.00  0.03  0.07  0.10   0.12   0.14   0.19 3658.44    1
#> mu0     13.07    0.09  4.46  2.72 10.66  13.56  16.01  20.82 2242.50    1
#> sigma0   7.33    0.06  3.97  1.48  4.48   6.69   9.48  16.56 3735.63    1
#> alpha    1.07    0.01  0.38  0.50  0.81   1.02   1.27   1.99 3474.06    1
```

## Nested example

Above, we have fitted the model only to a subset of trials in one
experimental condition. What is novel in `RStanTVA` is the possibility
to specify all kinds of hierarchical and nested parameter structures. At
first, we are now creating a subset of *all* trials of subject 20, which
includes trials in both conditions `high` and `low`:

``` r
tva_data_nested <- tva_recovery %>% filter(subject == 20)
tva_data_nested
#> # A tibble: 234 × 9
#> # Groups:   subject [1]
#>    subject trial     T condition    nD S[,1]  [,2] D[,1] R[,1] true_values$t[,1]
#>      <int> <int> <dbl> <fct>     <int> <int> <int> <int> <int>             <dbl>
#>  1      20     1  16.6 low           0     1     1     0     0              39.3
#>  2      20     2  33.3 high          0     1     1     0     1              34.5
#>  3      20     3  16.6 high          0     1     1     0     0              41.4
#>  4      20     4 150   high          4     1     1     1     0              56.8
#>  5      20     5  16.6 low           0     1     1     0     0              99.0
#>  6      20     6 100   low           0     1     1     0     1              29.0
#>  7      20     7 150   low           0     1     1     0     0              72.7
#>  8      20     8 150   low           0     1     1     0     1              28.6
#>  9      20     9 150   low           2     1     1     1     0              33.7
#> 10      20    10  16.6 low           0     1     1     0     0              79.0
#> # ℹ 224 more rows
#> # ℹ 13 more variables: S[3:6] <int>, D[2:6] <int>, R[2:6] <int>,
#> #   true_values$t[2:6] <dbl>, true_values$v <dbl[,6]>, $t0 <dbl>, $K <int>,
#> #   $C <dbl>, $alpha <dbl>, $w <dbl[,6]>, $mu0 <dbl>, $sigma0 <dbl>,
#> #   $pK <dbl[,7]>
```

We can now define a model for partial report of 6 display locations with
Gaussian $t_0$ and a free $K$ distribution, where we let $C$ vary
between conditions:

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

The “regression” formula for $C$ above will add another “layer” to the
model, which implements $C$ in trial $i$, $C_i$, as:

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
 *  This is a StanTVA program, generated with RStanTVA (v0.2.0) Please cite as:
 *  
 *  Rabe M, Kyllingsbæk S (2025). _RStanTVA: TVA Models in 'Stan' using 'R'
 *  and 'StanTVA'_. R package version 0.2.0,
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

You can do the same with any other model parameter and with the same
flexibility that you may be used to from `lme4` or `brms`. This also
includes the possibility of continuous covariates to directly model
effects of all kinds of biological, physical, or psychological metrics.
If you do not define a such a formula for any given parameter, it will
be fitted as a single invariant parameter across all trials. This can be
useful if you want to keep other parameters strictly constant across
conditions, e.g., in order to borrow statistical power.

We can fit `tva_model_nested` to the `tva_data_nested` using
maximum-likelihood estimation (MLE) or Bayesian parameter sampling in
the exactly same way as above:

``` r
tva_fit_nested <- sampling(tva_model_nested, tva_data_nested)
```

``` r
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
#>         mean se_mean   sd 2.5%  25%   50%   75% 97.5%   n_eff Rhat
#> w[1]    0.17    0.00 0.02 0.13 0.15  0.17  0.18  0.20 6389.54    1
#> w[2]    0.21    0.00 0.02 0.17 0.20  0.21  0.23  0.26 5931.22    1
#> w[3]    0.13    0.00 0.02 0.11 0.12  0.13  0.14  0.17 6396.66    1
#> w[4]    0.14    0.00 0.02 0.11 0.12  0.13  0.15  0.17 6375.33    1
#> w[5]    0.17    0.00 0.02 0.14 0.16  0.17  0.18  0.21 6887.20    1
#> w[6]    0.18    0.00 0.02 0.14 0.17  0.18  0.19  0.22 7444.45    1
#> pK[1]   0.10    0.00 0.02 0.06 0.08  0.10  0.11  0.14 6850.11    1
#> pK[2]   0.10    0.00 0.02 0.06 0.09  0.10  0.11  0.14 6559.53    1
#> pK[3]   0.21    0.00 0.03 0.15 0.19  0.21  0.23  0.28 5605.04    1
#> pK[4]   0.27    0.00 0.04 0.20 0.24  0.26  0.29  0.34 6691.55    1
#> pK[5]   0.11    0.00 0.03 0.07 0.09  0.11  0.13  0.16 5173.85    1
#> pK[6]   0.12    0.00 0.03 0.07 0.10  0.12  0.14  0.18 5815.02    1
#> pK[7]   0.10    0.00 0.03 0.05 0.08  0.10  0.11  0.16 4411.66    1
#> mu0    10.71    0.07 3.90 1.32 8.70 11.25 13.30 16.91 2789.80    1
#> sigma0  6.27    0.05 3.46 1.03 3.81  5.82  8.22 14.30 4072.46    1
#> alpha   1.01    0.00 0.28 0.55 0.83  0.98  1.17  1.63 5856.68    1
#> Warning: Unknown or uninitialised column: `param`.
#> 
#> Fixed effects:
#>                  mean se_mean   sd  2.5%  25%   50%  75% 97.5%   n_eff Rhat
#> C_Intercept      4.60    0.01 0.29  4.06  4.4  4.59 4.78  5.22 2959.86    1
#> C_conditionhigh -0.03    0.00 0.27 -0.56 -0.2 -0.02 0.15  0.50 4063.37    1
```

## Hierarchical example

If you have the necessary computational ressources, you can also fit a
single hierarchical model to the entire data set, where you can let the
parameters vary not only between conditions but even between, e.g.,
subjects, experimenters, locations, devices, etc. As for the previous
model, you may also choose to only let specific parameters vary, or just
all of them. You can specify combine model covariance matrices to
capture correlations between parameters (as is done below for $C$ and
$alpha$):

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

The power of hierarchical modelling is especially apparent when dealing
with sparse data. Therefore, we are here only looking at the first 100
trials of the first 10 subjects:

``` r
tva_hierarchical_subset <- tva_recovery %>% filter(subject <= 10 & trial <= 100)
```

``` r
tva_hierarchical_fit <- sampling(tva_hierarchical_model, tva_hierarchical_subset)
```

``` r
tva_hierarchical_fit
#> StanTVA model with 6 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
#> 
#> Model configuration:
#> formula = list(log(C) ~ 1 + condition + (1 + condition | C_alpha | subject), log(alpha) ~ 1 + (1 | C_alpha | subject), log(pK) ~ 1 + (1 | subject), mu0 ~ 1 + (1 | subject), log(sigma0) ~ 1 + (1 | subject), log(w) ~ 1 + (1 | subject))
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
#> priors = prior(normal(0, 0.07), dpar = C) + prior(normal(4, 0.2), coef = Intercept, dpar = C) + prior(normal(0, 0.07), dpar = alpha) + prior(normal(-0.2, 0.1), coef = Intercept, dpar = alpha) + prior(normal(0, 0.03), dpar = pK) + prior(normal(0, 0.1), coef = Intercept, dpar = pK) + prior(normal(0, 5), dpar = mu0) + prior(normal(30, 15), coef = Intercept, dpar = mu0) + prior(normal(0, 0.04), dpar = sigma0) + prior(normal(0, 0.2), coef = Intercept, dpar = sigma0) + prior(normal(0, 0.05), dpar = w) + prior(normal(0,      0.1), coef = Intercept, dpar = w) + prior(normal(0, 10), class = sd, coef = Intercept, dpar = mu0) + prior(normal(0, 5), class = sd, dpar = mu0)
#> sanity_checks = TRUE
#> debug_neginf_loglik = FALSE
#> 
#> Fixed effects:
#>                   mean se_mean   sd  2.5%   25%   50%   75% 97.5%   n_eff Rhat
#> C_Intercept       4.25    0.00 0.08  4.09  4.19  4.24  4.30  4.41 3634.74    1
#> C_conditionhigh   0.01    0.00 0.06 -0.10 -0.03  0.01  0.05  0.13 3818.13    1
#> alpha_Intercept  -0.25    0.00 0.08 -0.42 -0.31 -0.25 -0.20 -0.10 4208.60    1
#> pK_Intercept[1]  -0.18    0.00 0.08 -0.35 -0.24 -0.18 -0.12 -0.02 3475.77    1
#> pK_Intercept[2]  -0.08    0.00 0.09 -0.25 -0.14 -0.08 -0.02  0.08 3117.08    1
#> pK_Intercept[3]   0.21    0.00 0.09  0.03  0.15  0.21  0.27  0.38 2297.50    1
#> pK_Intercept[4]   0.11    0.00 0.09 -0.06  0.05  0.11  0.17  0.28 5884.16    1
#> pK_Intercept[5]  -0.03    0.00 0.09 -0.21 -0.08 -0.03  0.03  0.14 4882.48    1
#> pK_Intercept[6]  -0.02    0.00 0.09 -0.21 -0.09 -0.02  0.04  0.16 6240.98    1
#> mu0_Intercept    10.87    0.09 3.01  5.09  8.87 10.90 12.80 16.78 1199.21    1
#> sigma0_Intercept  0.02    0.00 0.20 -0.37 -0.10  0.02  0.16  0.40 4652.04    1
#> w_Intercept[1]   -0.06    0.00 0.06 -0.19 -0.11 -0.06 -0.02  0.07 4082.87    1
#> w_Intercept[2]    0.04    0.00 0.06 -0.08 -0.01  0.04  0.08  0.16 3846.57    1
#> w_Intercept[3]    0.01    0.00 0.06 -0.11 -0.03  0.01  0.06  0.14 4923.33    1
#> w_Intercept[4]    0.06    0.00 0.07 -0.08  0.02  0.06  0.10  0.18 1321.66    1
#> w_Intercept[5]    0.03    0.00 0.06 -0.08 -0.01  0.03  0.08  0.16 4374.83    1
#> 
#> Hyperparameters on random effects (subject level, N = 10):
#>                                                       mean se_mean   sd  2.5%
#> sd(C_subject_Intercept)                               0.07    0.00 0.05  0.00
#> sd(C_subject_conditionhigh)                           0.04    0.00 0.03  0.00
#> sd(alpha_subject_Intercept)                           0.08    0.00 0.06  0.00
#> sd(pK_subject_Intercept[1])                           0.13    0.01 0.09  0.01
#> sd(pK_subject_Intercept[2])                           0.08    0.01 0.06  0.00
#> sd(pK_subject_Intercept[3])                           0.13    0.01 0.09  0.01
#> sd(pK_subject_Intercept[4])                           0.10    0.00 0.07  0.01
#> sd(pK_subject_Intercept[5])                           0.08    0.00 0.06  0.01
#> sd(pK_subject_Intercept[6])                           0.08    0.00 0.06  0.00
#> sd(mu0_subject_Intercept)                             8.52    0.06 2.64  4.50
#> sd(sigma0_subject_Intercept)                          0.08    0.01 0.06  0.00
#> sd(w_subject_Intercept[1])                            0.08    0.00 0.06  0.01
#> sd(w_subject_Intercept[2])                            0.12    0.00 0.07  0.01
#> sd(w_subject_Intercept[3])                            0.08    0.00 0.05  0.01
#> sd(w_subject_Intercept[4])                            0.12    0.00 0.07  0.01
#> sd(w_subject_Intercept[5])                            0.08    0.00 0.05  0.01
#> cor(C_subject_Intercept,C_subject_conditionhigh)     -0.01    0.01 0.34 -0.65
#> cor(C_subject_Intercept,alpha_subject_Intercept)      0.00    0.01 0.35 -0.64
#> cor(C_subject_conditionhigh,alpha_subject_Intercept)  0.00    0.01 0.35 -0.67
#>                                                        25%   50%   75% 97.5%
#> sd(C_subject_Intercept)                               0.03  0.07  0.11  0.20
#> sd(C_subject_conditionhigh)                           0.02  0.03  0.06  0.11
#> sd(alpha_subject_Intercept)                           0.03  0.07  0.11  0.22
#> sd(pK_subject_Intercept[1])                           0.06  0.12  0.19  0.31
#> sd(pK_subject_Intercept[2])                           0.03  0.06  0.11  0.23
#> sd(pK_subject_Intercept[3])                           0.06  0.12  0.20  0.33
#> sd(pK_subject_Intercept[4])                           0.05  0.09  0.14  0.25
#> sd(pK_subject_Intercept[5])                           0.03  0.06  0.11  0.23
#> sd(pK_subject_Intercept[6])                           0.04  0.07  0.12  0.21
#> sd(mu0_subject_Intercept)                             6.63  8.13 10.03 14.47
#> sd(sigma0_subject_Intercept)                          0.03  0.07  0.12  0.22
#> sd(w_subject_Intercept[1])                            0.03  0.07  0.12  0.21
#> sd(w_subject_Intercept[2])                            0.06  0.11  0.16  0.26
#> sd(w_subject_Intercept[3])                            0.04  0.07  0.11  0.19
#> sd(w_subject_Intercept[4])                            0.06  0.12  0.17  0.26
#> sd(w_subject_Intercept[5])                            0.04  0.07  0.11  0.20
#> cor(C_subject_Intercept,C_subject_conditionhigh)     -0.26 -0.01  0.25  0.65
#> cor(C_subject_Intercept,alpha_subject_Intercept)     -0.25 -0.01  0.25  0.67
#> cor(C_subject_conditionhigh,alpha_subject_Intercept) -0.26  0.00  0.25  0.66
#>                                                        n_eff Rhat
#> sd(C_subject_Intercept)                               243.44 1.02
#> sd(C_subject_conditionhigh)                           217.85 1.03
#> sd(alpha_subject_Intercept)                           212.46 1.02
#> sd(pK_subject_Intercept[1])                           135.83 1.03
#> sd(pK_subject_Intercept[2])                           114.19 1.02
#> sd(pK_subject_Intercept[3])                           167.43 1.03
#> sd(pK_subject_Intercept[4])                           204.51 1.01
#> sd(pK_subject_Intercept[5])                           180.01 1.00
#> sd(pK_subject_Intercept[6])                           221.09 1.01
#> sd(mu0_subject_Intercept)                            2242.84 1.00
#> sd(sigma0_subject_Intercept)                          135.26 1.03
#> sd(w_subject_Intercept[1])                            174.37 1.01
#> sd(w_subject_Intercept[2])                            240.97 1.01
#> sd(w_subject_Intercept[3])                            249.55 1.01
#> sd(w_subject_Intercept[4])                            204.04 1.01
#> sd(w_subject_Intercept[5])                            300.78 1.01
#> cor(C_subject_Intercept,C_subject_conditionhigh)     3201.99 1.00
#> cor(C_subject_Intercept,alpha_subject_Intercept)     1814.57 1.00
#> cor(C_subject_conditionhigh,alpha_subject_Intercept) 3113.59 1.00
```

## Parallelization

There are two ways to speed up the parameter fitting: (1) Parallel
running of HMC chains, and (2) within-chain parallelization. The
technical details are described in more detail in the documenation of
Stan and `rstan`. You can use both, either, or none.

### Parallel chains (only relevant for `sampling()`)

By default, `rstan` will fit 4 chains sequentially. You can choose to
have them run in parallel instead by specifying the argument `cores` to
the number of chains to run in parallel. For example, to run all 4
chains in parallel, you can use:

``` r
fit <- sampling(tva_model, tva_data, chains = 4, cores = 4)
```

This may not work well in some cases, e.g., when run during knitting an
Rmd script.

### Within-chain parallelization (relevant for both `optimizing()` and `sampling()`)

You can also have the actual likelihood computation run in parallel,
which will calculate the likelihood for all trials independently,
scattered across separate *threads*, and ultimately aggregate them. We
are using Stan’s `map_rect()` function for this, which can even be used
for workload distribution across several compute nodes using MPI.

To take advantage of this, you need to explicitly specify the number of
CPUs **before** compiling the model. This is because
`stantva_model()`/`stantva_code()` will need to produce slightly
different (and more complex) Stan code when using this feature. For
example, if you want to use 8 CPUs:

``` r
rstan_options(threads_per_chain = 8)
tva_model <- stantva_model(...)
fit <- sampling(tva_model, tva_data)
```

If you are unsure how many CPUs are available on the target machine, you
can use `parallel::detectCores()` to determine that number:

``` r
rstan_options(threads_per_chain = parallel::detectCores())
```

Note that, for very complex models, you may notice average CPU loads
significantly below 100%. This is because some computational steps in
Stan cannot fully take advantage of the parallelization.

### Combining parallel chains and within-chain parallelization

As mentioned above, both can be combined but you should keep in mind
that the total number of requested CPUs is given by
$\textrm{threads}\times\textrm{chains}$ and may freeze your machine if
the average total CPU load exceeds the available ressources. In the case
of very complex models, where within-chain parallelization will hardly
achieve 100% on-average CPU loads, this may still be a good idea.

If you do not want to risk to overwhelm your machine, the number of
threads should always lie between the number of total CPUs
(`parallel::detectCores()`), and that number divided by the number of
chains, depending on the complexity of the model:

``` r
total_chains <- 4
parallel_chains <- parallel_chains # can also be lower
nthreads_per_chain <- ceiling(parallel::detectCores()/parallel_chains) # for simpler models
nthreads_per_chain <- parallel::detectCores() # for very complex models
rstan_options(threads_per_chain = nthreads_per_chain)
tva_model <- stantva_model(...)
fit <- sampling(tva_model, tva_data, chains = total_chains, cores = parallel_chains)
```
