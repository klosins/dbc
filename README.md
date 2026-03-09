
<!-- badges: start -->

[![R-CMD-check](https://github.com/klosins/dbc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/klosins/dbc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# dbc

## Dynamic bias correction for panel data estimators

This package implements the dynamic biases correction estimator as per
Klosin, S, Dynamic Biases of Static Panel Data Estimators (forthcoming).
This allows the estimation of treatment effects in panel fixed effects
models where the outcome is dynamic.

## Installation

``` r
# Install the CRAN version
pak::pak("dbc")

# Or install the latest version on github
pak::pak("github::klosins/dbc")
```

## Usage

``` r
library(dbc)
```

The follows shows example usage of the DBC under three general types of
model.

The first is the **exogenous** treatment model, where the outcome is:
$$Y_{it} = \alpha_{i} + \rho_{1} Y_{t-1} + \tau D_{t} + \beta_{1} X_{1} +\varepsilon_{it}$$

``` r


# Exogenous treatment (no lag_y in treatment equation)
data_exog <- DGP(N = 500, N_T = 4, rho1 = 0.2, rho2 = 0, tau = 0.5,n_X1 = 1, n_X2 = 0)


fit_exog <- dbc(
    outcome_fml = y ~ lag_y + D,
    lag_y       = "lag_y",
    treatment   = "D",
    panel_id    = "panel_id",
    time_id     = "time",
    data        = data_exog
)
#> Warning in dbc(outcome_fml = y ~ lag_y + D, lag_y = "lag_y", treatment = "D", : Estimated |phi| = 1.3463 >= 1. Stationarity assumption (Assumption 5, eq 35) violated. Bias correction may be
#> unreliable.


summary(fit_exog)
#> Dynamic Bias-Corrected Estimator
#> 
#> Call:
#> dbc(outcome_fml = y ~ lag_y + D, lag_y = "lag_y", treatment = "D",     panel_id = "panel_id", time_id = "time", data = data_exog)
#> 
#> Outcome equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y  1.34626    0.02326  57.883   <2e-16 ***
#> D      0.85785    0.08741   9.814   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> N = 500  T = 3  phi = 1.3463  GMM: (converged)
```

The second, is the **endogenous** treatment model, where the treatment
effect is also a function of the previous values of the outcome:
$$Y_{it} = \alpha_{i} + \rho_{1} Y_{t-1} + \tau D_{t} \beta_{1}X_{1}  +\varepsilon_{it}$$
$$D_{it} = c_{i} + \rho_{2} Y_{t-1}  + \beta_{2} X_{2} + u_{it}$$

``` r

set.seed(42)

# Endogenous treatment (lag_y in treatment equation)
data_endo <- DGP(N = 500, N_T = 4, rho1 = 0.2, rho2 = 0.3, tau = 0.5,n_X1 = 1, n_X2 = 1)

fit_endo <- dbc(
    outcome_fml   = y ~ lag_y + D,
    treatment_fml = D ~ lag_y,
    lag_y         = "lag_y",
    treatment     = "D",
    panel_id      = "panel_id",
    time_id       = "time",
    data          = data_endo
)

summary(fit_endo)
#> Dynamic Bias-Corrected Estimator
#> 
#> Call:
#> dbc(outcome_fml = y ~ lag_y + D, treatment_fml = D ~ lag_y, lag_y = "lag_y",     treatment = "D", panel_id = "panel_id", time_id = "time",     data = data_endo)
#> 
#> Outcome equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y  0.19012    0.02422    7.85 4.15e-15 ***
#> D      0.48392    0.03416   14.17  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Treatment equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y  0.32723    0.02197   14.89   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> N = 500  T = 3  phi = 0.3485  GMM: (converged)
```

Last, the model that allows for interactions between the treatment
variable and exogenous variables:
$$Y_{it} = \alpha_{i} + \rho_{1} Y_{t-1} + \tau_{1} W_{1} \times D_{t} + \tau_{2} W_{2} \times D_{t} \beta_{1}X_{1}  +\varepsilon_{it}$$
$$D_{it} = c_{i} + \rho_{2} Y_{t-1}  + \beta_{2} X_{2} + u_{it}$$

``` r


# With moderators/interactions
data_interact <- DGP(
    N = 500, N_T = 4,
    rho1 = 0.2, rho2 = 0.3,
    tau = c(0.8, 1.2),
    n_X1 = 1, n_X2 = 1, n_W = 2
)


fit_interact <- dbc(
    outcome_fml   = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
    treatment_fml = D ~ lag_y + X2_1,
    lag_y         = "lag_y",
    treatment     = "D",
    panel_id      = "panel_id",
    time_id       = "time",
    data          = data_interact
)

summary(fit_interact)
#> Dynamic Bias-Corrected Estimator
#> 
#> Call:
#> dbc(outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1, treatment_fml = D ~     lag_y + X2_1, lag_y = "lag_y", treatment = "D", panel_id = "panel_id",     time_id = "time", data = data_interact)
#> 
#> Outcome equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y 0.197031   0.003657   53.88   <2e-16 ***
#> D:W_1 0.793227   0.004728  167.77   <2e-16 ***
#> D:W_2 1.199933   0.004571  262.50   <2e-16 ***
#> X1_1  0.992106   0.032193   30.82   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Treatment equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y 0.297533   0.003347   88.90   <2e-16 ***
#> X2_1  0.998390   0.029439   33.91   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> N = 500  T = 3  phi = 0.2007  GMM: (converged)
```

### S3 methods

Standard model sumamry functions on the dbc output are provided.

Returning the bias-corrected coefficients:

``` r
coef(fit_endo)
#>          lag_y              D treat_eq:lag_y 
#>      0.1901217      0.4839166      0.3272342
```

Returning the variance-convariance matrix on the bias corrected
coefficients:

``` r
vcov(fit_endo)
#>                        lag_y             D treat_eq:lag_y
#> lag_y           5.865491e-04 -0.0003524152  -4.479349e-05
#> D              -3.524152e-04  0.0011668063   2.971201e-04
#> treat_eq:lag_y -4.479349e-05  0.0002971201   4.827017e-04
```

Standard confidence intervals:

``` r
confint(fit_endo, level = 0.90)
#>                       5%       95%
#> lag_y          0.1502853 0.2299580
#> D              0.4277308 0.5401024
#> treat_eq:lag_y 0.2910960 0.3633724
```

The number of observations:

``` r
nobs(fit_endo)
#> [1] 1500
```

### broom / modelsummary

The package has methods to interface with the `broom` and `modelsummary`
packages.

Create a tidy version of the bias-corrected coefficients:

``` r
library(broom)

tidy(fit_endo, conf.int = TRUE, conf.level = 0.95)
#>       group  term  estimate  std.error statistic      p.value  conf.low conf.high
#> 1   Outcome lag_y 0.1901217 0.02421878  7.850176 4.154532e-15 0.1426537 0.2375896
#> 2   Outcome     D 0.4839166 0.03415855 14.166779 1.470947e-45 0.4169670 0.5508661
#> 3 Treatment lag_y 0.3272342 0.02197047 14.894272 3.590735e-50 0.2841729 0.3702955
```

`glance` shows other model information, like number of observations.

``` r
glance(fit_endo)
#>         N N_T nobs   phi converged
#> lag_y 500   3 1500 0.348      TRUE
```

We can produce tables from `modelsummary`. Using the
`shape = group + term ~ model` argument groups the table by whether the
coefficients are in the outcome or treatment equation.

``` r
library(modelsummary)

# Multiple models side-by-side
modelsummary(
    list(Exogenous = fit_exog, Endogenous = fit_endo),
    shape = group + term ~ model
)
```

|           |          | Exogenous | Endogenous |
|-----------|----------|-----------|------------|
| Outcome   | lag_y    | 1.346     | 0.190      |
|           |          | (0.023)   | (0.024)    |
|           | D        | 0.858     | 0.484      |
|           |          | (0.087)   | (0.034)    |
| Treatment | lag_y    |           | 0.327      |
|           |          |           | (0.022)    |
|           | Num.Obs. | 1500      | 1500       |
|           | N        | 500       | 500        |
|           | N_T      | 3         | 3          |
|           | phi      | 1.35      | 0.348      |
