
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RStanTVA

<!-- badges: start -->

[![R](https://github.com/mmrabe/RStanTVA/actions/workflows/rpkg.yml/badge.svg)](https://github.com/mmrabe/RStanTVA/actions/workflows/rpkg.yml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/RStanTVA)](https://CRAN.R-project.org/package=RStanTVA)
<!-- badges: end -->

RStanTVA is an R package containing the StanTVA library and numerous
convenience functions for generating, compiling, fitting, and analyzing
(Stan)TVA models.

## Installation

The latest version of the package can be installed from CRAN using:

``` r
install.packages("RStanTVA")
```

Alternatively, you can install the development version of RStanTVA from
[GitHub](https://github.com/mmrabe/RStanTVA) with:

``` r
remotes::install_github("mmrabe/RStanTVA")
```

## Loading the library

After installing, RStanTVA can be loaded like any other R package. For
data wrangling, you may also find it helpful to load dplyr and tidyr (or
the complete tidyverse), which are also used in this readme.

``` r
library(RStanTVA)
library(tidyverse)
```

Loading RStanTVA will also load rstan and its dependencies, which may
produce additional messages in the console. Regarding rstan’s note on
multi-threading, please see the Parallelization section at the end of
this readme.

## Example data set

The functions in this package use a relatively simplistic data format.
For an experiment with $n$ trials and $m$ display locations, the data
should be represented as a `list` or `tibble` with the following
columns/elements:

- `S`: An $n \times m$ matrix with values `TRUE` or `1` wherever a
  stimulus was displayed, otherwise `FALSE` or `0`. If stimuli were
  displayed in all $m$ display locations in all $n$ trials, this would
  simply be `matrix(TRUE,n,m)`.
- `R`: An $n \times m$ matrix with values `TRUE` or `1` for every
  correctly reported *target* location, otherwise `FALSE` or `0`. This
  must be `FALSE`/`0` where `S` is `FALSE`/`0` or `D` (see below) is
  `TRUE`/`1`, otherwise will throw an error.
- `D`: An $n \times m$ matrix with values `TRUE` or `1` for every
  location, otherwise `FALSE` or `0`. This must be `FALSE`/`0` where `S`
  is `FALSE`/`0` or `R` (see below) is `TRUE`/`1`, otherwise will throw
  an error. `D` is optional if *all* trials are whole report.
  Whole-report trials can also be included by setting all `D` of the
  respective trials to `FALSE`/`0`, since a partial report without
  distractors is effectively a whole report.
- `T`: An $n$-length numeric vector of exposure durations (in
  milliseconds).

``` r
data("tva_recovery")
data("tva_recovery_true_params")
```

`RStanTVA` comes with a CombiTVA data set `tva_recovery` of 50 simulated
subjects with 234 trials each. For each subject and trial, the true
underlying parameters are known exactly, which makes it very useful for
demonstrating the functionality of this package:

Of the 234 trials, half were carried out under the `high` condition, and
the other half under the `low` condition. Imagine that those are maybe
different levels of luminance, some sort of cognitive load or so. The
concrete manipulation is not really relevant but we *do* know that on
average, $C_\mathrm{high}=100.0>C_\mathrm{low}=80.0$, and that across
subjects, $\mathrm{cor}\left(\log C,\log\alpha\right)=-0.3$.

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
#>  1      20     1  16.6 high          0     1     1     0     0              37.0
#>  2      20     5 100   high          0     1     1     0     0              40.5
#>  3      20     7 150   high          2     1     1     1     0             302. 
#>  4      20    11 150   high          3     1     1     0     0             123. 
#>  5      20    12  50   high          0     1     1     0     0              36.8
#>  6      20    15  33.3 high          0     1     1     0     0              75.3
#>  7      20    16 100   high          0     1     1     0     0             115. 
#>  8      20    17  16.6 high          0     1     1     0     0              74.6
#>  9      20    18 150   high          2     1     1     1     0              71.3
#> 10      20    19 100   high          0     1     1     0     1              49.9
#> # ℹ 107 more rows
#> # ℹ 14 more variables: S[3:6] <int>, D[2:6] <int>, R[2:6] <int>,
#> #   true_values$t[2:6] <dbl>, true_values$v <dbl[,6]>, $t0 <dbl>, $K <dbl>,
#> #   $C <dbl>, $alpha <dbl>, $w <dbl[,2]>, $mu0 <dbl>, $sigma0 <dbl>, $mK <dbl>,
#> #   $sK <dbl>
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
  w_mode = "regions",
  regions = list(1:3, 4:6),
  t0_mode = "gaussian",
  K_mode = "probit"
)
```

If we wanted to look at the generated Stan code, this would be it:

``` stan
/************************************************************************************
 *  StanTVA
 *  =======
 *  This is a StanTVA program, generated with RStanTVA (v0.3.0) Please cite as:
 *  
 *  Rabe M, Kyllingsbæk S (2025). _RStanTVA: TVA Models in 'Stan' using 'R'
 *  and 'StanTVA'_. R package version 0.3.0,
 *  <https://github.com/mmrabe/RStanTVA>.
 *  
 *  Configuration
 *  =============
 *  - formula = NULL
 *  - locations = 6
 *  - task = pr
 *  - regions = list(1:3, 4:6)
 *  - C_mode = equal
 *  - w_mode = regions
 *  - t0_mode = gaussian
 *  - K_mode = probit
 *  - max_K = 6
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
    #include probitK.stan
    #include gaussiant0.stan
    vector calculate_v(data int nS, data array[] int S, data array[] int D, vector r, real C, real alpha) {
        vector[6] w;
        w[1] = r[1]/3.0;
        w[2] = r[1]/3.0;
        w[3] = r[1]/3.0;
        w[4] = r[2]/3.0;
        w[5] = r[2]/3.0;
        w[6] = r[2]/3.0;
        vector[6] s = rep_vector(C, 6);
        array[nS] int Ss = get_matches(S);
        vector[6] w_alpha = w;
        for(i in 1:6) if(D[i]) w_alpha[i] *= alpha;
        vector[nS] v = s[Ss] .* w_alpha[Ss] / sum(w_alpha[Ss]);
        for(i in 1:nS) if(v[i] < machine_precision()) v[i] = machine_precision();
        return v/1000.0;
    }
    real log_lik_single(data array[] int S, data array[] int D, data array[] int R, data int nS, data real T, vector r, real C, real alpha, real mK, real sK, real mu0, real sigma0) {
        real log_lik;
        vector[nS] v = calculate_v(nS, S, D, to_vector(r), C, alpha);
        log_lik = tva_pr_log(R, S, D, T, [mu0, sigma0]', [mK, sK]', v);
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
    int max_K;
    max_K = 6;
    array[N] int nS;
    for(i in 1:N) nS[i] = sum(S[i,]);
    int total_nS = sum(nS);
}
parameters {
    real<lower=machine_precision()> C;
    simplex[2] r;
    real mK;
    real<lower=machine_precision()> sK;
    real mu0;
    real<lower=machine_precision()> sigma0;
    real<lower=machine_precision()> alpha;
}
model {
    C ~ lognormal(4.5, 0.6);
    r[:1]/r[2] ~ lognormal(0, 0.2);
    mK ~ normal(3.5, 0.5);
    sK ~ lognormal(-0.5, 0.25);
    mu0 ~ normal(20, 15);
    sigma0 ~ lognormal(2, 0.6);
    alpha ~ lognormal(-0.4, 0.6);
    // likelihood (only if prior != 0)
    if(target() != negative_infinity()) {
        for(i in 1:N) target += log_lik_single(S[i], D[i], R[i], nS[i], T[i], to_vector(r), C, alpha, mK, sK, mu0, sigma0);
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
tva_fit_mle$par[c("C","alpha","mu0","sigma0","mK","sK")]
#>          C      alpha        mu0     sigma0         mK         sK 
#> 95.8549318  0.6201870 21.4276180  7.7779863  3.0318789  0.4844716
```

… or using Bayesian HMC sampling:

``` r
tva_fit <- sampling(tva_model, tva_data)
```

``` r
tva_fit
#> StanTVA model with 7 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
#> 
#> Model configuration:
#> locations = 6
#> task = "pr"
#> regions = list(1:3, 4:6)
#> C_mode = "equal"
#> w_mode = "regions"
#> t0_mode = "gaussian"
#> K_mode = "probit"
#> max_K = 6
#> parallel = TRUE
#> save_log_lik = FALSE
#> sanity_checks = TRUE
#> debug_neginf_loglik = FALSE
#> Warning: Unknown or uninitialised column: `param`.
#> 
#> Global parameters:
#>          mean se_mean    sd  2.5%   25%    50%    75%  97.5%   n_eff Rhat
#> C      111.38    0.85 36.19 63.00 86.23 104.45 129.01 200.13 1814.83    1
#> r[1]     0.49    0.00  0.03  0.43  0.47   0.49   0.51   0.55 4843.93    1
#> r[2]     0.51    0.00  0.03  0.45  0.49   0.51   0.53   0.57 4843.93    1
#> mK       3.04    0.00  0.09  2.86  2.97   3.04   3.10   3.21 3038.93    1
#> sK       0.51    0.00  0.07  0.38  0.46   0.50   0.55   0.65 4404.91    1
#> mu0     22.49    0.10  4.62 13.46 19.44  22.50  25.58  31.40 2036.00    1
#> sigma0  10.19    0.07  3.85  3.73  7.43   9.91  12.45  18.90 3483.68    1
#> alpha    0.67    0.00  0.17  0.39  0.55   0.65   0.77   1.05 4751.24    1
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
#>  1      20     1  16.6 high          0     1     1     0     0              37.0
#>  2      20     2 150   low           4     1     1     1     0             298. 
#>  3      20     3 150   low           0     1     1     0     0             126. 
#>  4      20     4 200   low           0     1     1     0     0              75.9
#>  5      20     5 100   high          0     1     1     0     0              40.5
#>  6      20     6 150   low           0     1     1     0     0              70.4
#>  7      20     7 150   high          2     1     1     1     0             302. 
#>  8      20     8 150   low           2     1     1     0     1             178. 
#>  9      20     9  50   low           0     1     1     0     1             146. 
#> 10      20    10  16.6 low           0     1     1     0     0              63.8
#> # ℹ 224 more rows
#> # ℹ 14 more variables: S[3:6] <int>, D[2:6] <int>, R[2:6] <int>,
#> #   true_values$t[2:6] <dbl>, true_values$v <dbl[,6]>, $t0 <dbl>, $K <dbl>,
#> #   $C <dbl>, $alpha <dbl>, $w <dbl[,2]>, $mu0 <dbl>, $sigma0 <dbl>, $mK <dbl>,
#> #   $sK <dbl>
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
  w_mode = "regions",
  regions = list(1:3, 4:6),
  t0_mode = "gaussian",
  K_mode = "probit"
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
 *  This is a StanTVA program, generated with RStanTVA (v0.3.0) Please cite as:
 *  
 *  Rabe M, Kyllingsbæk S (2025). _RStanTVA: TVA Models in 'Stan' using 'R'
 *  and 'StanTVA'_. R package version 0.3.0,
 *  <https://github.com/mmrabe/RStanTVA>.
 *  
 *  Configuration
 *  =============
 *  - formula = list(log(C) ~ 1 + condition)
 *  - locations = 6
 *  - task = pr
 *  - regions = list(1:3, 4:6)
 *  - C_mode = equal
 *  - w_mode = regions
 *  - t0_mode = gaussian
 *  - K_mode = probit
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
    #include probitK.stan
    #include gaussiant0.stan
    vector calculate_v(data int nS, data array[] int S, data array[] int D, vector r, real alpha, real C) {
        vector[6] w;
        w[1] = r[1]/3.0;
        w[2] = r[1]/3.0;
        w[3] = r[1]/3.0;
        w[4] = r[2]/3.0;
        w[5] = r[2]/3.0;
        w[6] = r[2]/3.0;
        vector[6] s = rep_vector(C, 6);
        array[nS] int Ss = get_matches(S);
        vector[6] w_alpha = w;
        for(i in 1:6) if(D[i]) w_alpha[i] *= alpha;
        vector[nS] v = s[Ss] .* w_alpha[Ss] / sum(w_alpha[Ss]);
        for(i in 1:nS) if(v[i] < machine_precision()) v[i] = machine_precision();
        return v/1000.0;
    }
    real log_lik_single(data array[] int S, data array[] int D, data array[] int R, data int nS, data real T, vector r, real alpha, real mK, real sK, real mu0, real sigma0, real C) {
        real log_lik;
        vector[nS] v = calculate_v(nS, S, D, to_vector(r), alpha, C);
        log_lik = tva_pr_log(R, S, D, T, [mu0, sigma0]', [mK, sK]', v);
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
    int max_K;
    max_K = 6;
    array[N] int nS;
    for(i in 1:N) nS[i] = sum(S[i,]);
    int total_nS = sum(nS);
    int M = M_C;
}
parameters {
    simplex[2] r;
    real<lower=machine_precision()> mK;
    real<lower=machine_precision()> sK;
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
    r[:1]/r[2] ~ lognormal(0, 0.2);
    mK ~ normal(3.5, 1);
    sK ~ lognormal(0, 1);
    mu0 ~ normal(20, 15);
    sigma0 ~ gamma(2, 0.2);
    alpha ~ lognormal(-0.4, 0.6);
    // likelihood (only if prior != 0)
    if(target() != negative_infinity()) {
        for(i in 1:N) target += log_lik_single(S[i], D[i], R[i], nS[i], T[i], to_vector(r), alpha, mK, sK, mu0, sigma0, C[i]);
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
#> StanTVA model with 7 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
#> 
#> Model configuration:
#> formula = list(log(C) ~ 1 + condition)
#> locations = 6
#> task = "pr"
#> regions = list(1:3, 4:6)
#> C_mode = "equal"
#> w_mode = "regions"
#> t0_mode = "gaussian"
#> K_mode = "probit"
#> max_K = 6
#> allow_guessing = FALSE
#> parallel = TRUE
#> save_log_lik = FALSE
#> sanity_checks = TRUE
#> debug_neginf_loglik = FALSE
#> 
#> Global parameters:
#>         mean se_mean   sd  2.5%   25%   50%   75% 97.5%   n_eff Rhat
#> r[1]    0.48    0.00 0.02  0.43  0.46  0.48  0.49  0.52 3818.53    1
#> r[2]    0.52    0.00 0.02  0.48  0.51  0.52  0.54  0.57 3818.53    1
#> mK      3.07    0.00 0.07  2.93  3.03  3.07  3.12  3.21 3170.68    1
#> sK      1.93    0.00 0.22  1.52  1.78  1.92  2.07  2.37 2928.26    1
#> mu0    19.25    0.07 3.43 13.00 16.79 19.14 21.55 26.43 2385.18    1
#> sigma0  7.63    0.07 3.46  1.62  5.17  7.38  9.82 15.26 2554.38    1
#> alpha   0.64    0.00 0.12  0.43  0.55  0.63  0.71  0.90 3612.55    1
#> 
#> Fixed effects:
#>                 mean se_mean   sd  2.5%   25%  50%  75% 97.5%   n_eff Rhat
#> C_Intercept     4.43    0.01 0.24  4.01  4.26 4.40 4.57  4.98 1907.77    1
#> C_conditionhigh 0.07    0.00 0.26 -0.44 -0.11 0.06 0.24  0.58 3062.82    1
```

## Hierarchical example

If you have the necessary computational ressources, you can also fit a
single hierarchical model to the entire data set, where you can let the
parameters vary not only between conditions but even between, e.g.,
subjects, experimenters, locations, devices, etc. As for the previous
model, you may also choose to only let specific parameters vary, or just
all of them. You can specify combine model covariance matrices to
capture correlations between parameters (as is done below for $C$ and
$\alpha$):

``` r


priors <-
  prior(normal(0,.07),dpar=C)+
  prior(normal(4,.2),dpar=C,coef=Intercept)+
  prior(normal(0,.07),dpar=alpha)+
  prior(normal(-0.2,.1),dpar=alpha,coef=Intercept)+
  prior(normal(3.5,.5),dpar=mK,coef=Intercept)+
  prior(normal(0,.1),dpar=mK)+
  prior(normal(-0.5,.5),dpar=sK,coef=Intercept)+
  prior(normal(0,.1),dpar=sK)+
  prior(normal(0,5),dpar=mu0)+
  prior(normal(30,15),dpar=mu0,coef=Intercept)+
  prior(normal(0,.04),dpar=sigma0)+
  prior(normal(0,.2),dpar=sigma0,coef=Intercept)+
  prior(normal(0,.05),dpar=r)+
  prior(normal(0,0.1),dpar=r,coef=Intercept)

tva_hierarchical_model <- stantva_model(
  formula = list(
    log(C) ~ 1 + condition + (1 + condition | C_alpha | subject),
    log(r) ~ 1 + (1 | subject),
    log(mK) ~ 1 + (1 | K | subject),
    log(sK) ~ 1 + (1 | K | subject),
    mu0 ~ 1 + (1 | subject),
    log(sigma0) ~ 1 + (1 | subject),
    log(alpha) ~ 1 + (1 | C_alpha | subject)
  ),
  locations = 6,
  task = "pr",
  w_mode = "regions",
  regions = list(1:3, 4:6),
  t0_mode = "gaussian",
  K_mode = "probit",
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
tva_hierarchical_fit
```

    #> StanTVA model with 7 free parameter(s), fitted with 4  chains, each with iter=2000; warmup=1000; thin=1
    #> 
    #> Model configuration:
    #> formula = list(log(C) ~ 1 + condition + (1 + condition | C_alpha | subject), log(alpha) ~ 1 + (1 | C_alpha | subject), mK ~ 1 + (1 | K | subject), log(sK) ~ 1 + (1 | K | subject), mu0 ~ 1 + (1 | subject), log(sigma0) ~ 1 + (1 | subject), log(r) ~ 1 + (1 | subject))
    #> locations = 6
    #> task = "pr"
    #> regions = list(1:3, 4:6)
    #> C_mode = "equal"
    #> w_mode = "regions"
    #> t0_mode = "gaussian"
    #> K_mode = "probit"
    #> max_K = 6
    #> parallel = TRUE
    #> save_log_lik = FALSE
    #> priors = prior(normal(0, 0.07), dpar = C) + prior(normal(4, 0.2), coef = Intercept, dpar = C) + prior(normal(0, 0.07), dpar = alpha) + prior(normal(-0.2, 0.1), coef = Intercept, dpar = alpha) + prior(normal(3.5, 0.5), coef = Intercept, dpar = mK) + prior(normal(0, 0.1), dpar = mK) + prior(normal(-0.5, 0.5), coef = Intercept, dpar = sK) + prior(normal(0, 0.1), dpar = sK) + prior(normal(0, 5), dpar = mu0) + prior(normal(30, 15), coef = Intercept, dpar = mu0) + prior(normal(0, 0.2), dpar = sigma0) + prior(normal(1.7,      0.3), coef = Intercept, dpar = sigma0) + prior(normal(0, 0.05), dpar = r) + prior(normal(0, 0.1), coef = Intercept, dpar = r) + prior(gamma(2, 8), class = sd) + prior(gamma(2, 5), class = sd, coef = Intercept) + prior(gamma(2, 0.5), class = sd, dpar = mu0) + prior(gamma(2, 0.15), class = sd, coef = Intercept, dpar = mu0)
    #> sanity_checks = TRUE
    #> debug_neginf_loglik = FALSE
    #> 
    #> Fixed effects:
    #>                   mean se_mean   sd  2.5%   25%   50%   75% 97.5%   n_eff Rhat
    #> C_Intercept       4.36    0.00 0.11  4.14  4.30  4.37  4.44  4.55 1393.47    1
    #> C_conditionhigh   0.04    0.00 0.06 -0.08  0.00  0.04  0.08  0.16 3736.56    1
    #> alpha_Intercept  -0.23    0.00 0.08 -0.39 -0.28 -0.23 -0.17 -0.06 2649.44    1
    #> mK_Intercept      3.42    0.01 0.30  2.83  3.22  3.42  3.61  4.03  762.12    1
    #> sK_Intercept     -0.79    0.00 0.10 -0.99 -0.84 -0.78 -0.72 -0.61 1120.87    1
    #> mu0_Intercept    19.58    0.08 2.74 14.12 17.88 19.56 21.34 25.10 1298.38    1
    #> sigma0_Intercept  1.64    0.00 0.22  1.17  1.50  1.65  1.80  2.05 2272.88    1
    #> r_Intercept[1]    0.02    0.00 0.07 -0.11 -0.02  0.02  0.07  0.15 2515.91    1
    #> 
    #> Hyperparameters on random effects (subject level, N = 10):
    #>                                                       mean se_mean   sd  2.5%
    #> sd(C_subject_Intercept)                               0.19    0.00 0.12  0.03
    #> sd(C_subject_conditionhigh)                           0.15    0.00 0.10  0.02
    #> sd(alpha_subject_Intercept)                           0.38    0.00 0.17  0.09
    #> sd(mK_subject_Intercept)                              1.21    0.00 0.23  0.85
    #> sd(sK_subject_Intercept)                              0.16    0.01 0.10  0.03
    #> sd(mu0_subject_Intercept)                             7.43    0.06 2.42  3.78
    #> sd(sigma0_subject_Intercept)                          0.31    0.01 0.21  0.05
    #> sd(r_subject_Intercept[1])                            0.22    0.00 0.09  0.08
    #> cor(C_subject_Intercept,C_subject_conditionhigh)      0.03    0.01 0.52 -0.89
    #> cor(C_subject_Intercept,alpha_subject_Intercept)      0.02    0.01 0.51 -0.88
    #> cor(C_subject_conditionhigh,alpha_subject_Intercept) -0.08    0.01 0.51 -0.90
    #> cor(mK_subject_Intercept,sK_subject_Intercept)       -0.09    0.01 0.54 -0.94
    #>                                                        25%   50%  75% 97.5%
    #> sd(C_subject_Intercept)                               0.11  0.17 0.25  0.48
    #> sd(C_subject_conditionhigh)                           0.08  0.13 0.21  0.41
    #> sd(alpha_subject_Intercept)                           0.26  0.37 0.49  0.76
    #> sd(mK_subject_Intercept)                              1.05  1.18 1.34  1.74
    #> sd(sK_subject_Intercept)                              0.09  0.14 0.21  0.41
    #> sd(mu0_subject_Intercept)                             5.72  7.02 8.75 13.11
    #> sd(sigma0_subject_Intercept)                          0.16  0.26 0.42  0.84
    #> sd(r_subject_Intercept[1])                            0.16  0.21 0.27  0.41
    #> cor(C_subject_Intercept,C_subject_conditionhigh)     -0.40  0.04 0.47  0.92
    #> cor(C_subject_Intercept,alpha_subject_Intercept)     -0.41  0.01 0.45  0.89
    #> cor(C_subject_conditionhigh,alpha_subject_Intercept) -0.49 -0.12 0.32  0.87
    #> cor(mK_subject_Intercept,sK_subject_Intercept)       -0.56 -0.11 0.35  0.88
    #>                                                        n_eff Rhat
    #> sd(C_subject_Intercept)                               642.83 1.00
    #> sd(C_subject_conditionhigh)                           589.38 1.00
    #> sd(alpha_subject_Intercept)                          1184.65 1.00
    #> sd(mK_subject_Intercept)                             3300.77 1.00
    #> sd(sK_subject_Intercept)                              356.65 1.01
    #> sd(mu0_subject_Intercept)                            1827.91 1.00
    #> sd(sigma0_subject_Intercept)                          231.15 1.01
    #> sd(r_subject_Intercept[1])                           1559.18 1.00
    #> cor(C_subject_Intercept,C_subject_conditionhigh)     1452.76 1.00
    #> cor(C_subject_Intercept,alpha_subject_Intercept)     1303.57 1.00
    #> cor(C_subject_conditionhigh,alpha_subject_Intercept) 1468.44 1.00
    #> cor(mK_subject_Intercept,sK_subject_Intercept)       1426.15 1.00

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
parallel_chains <- total_chains # can also be lower
nthreads_per_chain <- ceiling(parallel::detectCores()/parallel_chains) # for simpler models
nthreads_per_chain <- parallel::detectCores() # for very complex models
rstan_options(threads_per_chain = nthreads_per_chain)
tva_model <- stantva_model(...)
fit <- sampling(tva_model, tva_data, chains = total_chains, cores = parallel_chains)
```
