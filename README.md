
<!-- badges: start -->

[![R-CMD-check](https://github.com/klosins/dbc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/klosins/dbc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# dbc

## Dynamic bias correction for panel data estimators

This package implements the dynamic biases correction estimator as per
Klosin, S, Dynamic Biases of Static Panel Data Estimators (forthcoming).
This allows the estimation of treatment effects in panel fixed effects
models where the outcome is dynamic. See below for usage and examples.

## Installation

``` r
# Install the latest version on github
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

and the treatment equation is:

$$D_{it} = u_{it}$$.

The following shows the simulation of this model using the `dbc::DGP`
data generating process function, and then the estimation of this model.

`summary(fit_exog)` shows the bias-corrected estimates. The reported
standard errors used underlying iid errors, and are not clustered.

The biased, OLS estimates of the model can be shown as
`fit_exog$biased_coefficients`. `str(fit_exog)` shows the full list of
objects returned with the `dbc` model.

``` r

set.seed(42)

# Exogenous treatment (no lag_y in treatment equation)
data_exog <- DGP(N = 500, N_T = 4, rho1 = 0.2, rho2 = 0, tau = 0.5, n_X1 = 1, n_X2 = 0)

fit_exog <- dbc(
    outcome_fml = y ~ lag_y + D + X1_1,
    lag_y       = "lag_y",
    treatment   = "D",
    panel_id    = "panel_id",
    time_id     = "time",
    data        = data_exog
)


summary(fit_exog)
#> Dynamic Bias-Corrected Estimator
#> 
#> Call:
#> dbc(outcome_fml = y ~ lag_y + D + X1_1, lag_y = "lag_y", treatment = "D",     panel_id = "panel_id", time_id = "time", data = data_exog)
#> 
#> Outcome equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y  0.21457    0.02129   10.08   <2e-16 ***
#> D      0.46353    0.03051   15.19   <2e-16 ***
#> X1_1   0.95184    0.03370   28.25   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> N = 500  T = 3  phi = 0.2146  GMM: (converged)
```

The second, is the **endogenous** treatment model, where the treatment
effect is also a function of the previous values of the outcome:

$$Y_{it} = \alpha_{i} + \rho_{1} Y_{t-1} + \tau D_{t} \beta_{1}X_{1}  +\varepsilon_{it}$$

$$D_{it} = c_{i} + \rho_{2} Y_{t-1}  + \beta_{2} X_{2} + u_{it}$$

``` r


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
#> lag_y  0.19601    0.02423    8.09 5.99e-16 ***
#> D      0.47204    0.03146   15.00  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Treatment equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y  0.28360    0.02233    12.7   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> N = 500  T = 3  phi = 0.3299  GMM: (converged)
```

Last, the **interaction** model, which allows for interactions between
the treatment variable and exogenous variables, $W$:

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
#> lag_y 0.201906   0.003954   51.07   <2e-16 ***
#> D:W_1 0.800101   0.004293  186.36   <2e-16 ***
#> D:W_2 1.201314   0.003830  313.65   <2e-16 ***
#> X1_1  1.017534   0.032734   31.09   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Treatment equation:
#>       Estimate Std. Error z value Pr(>|z|)    
#> lag_y  0.30295    0.00296  102.35   <2e-16 ***
#> X2_1   1.00565    0.03086   32.58   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> N = 500  T = 3  phi = 0.2051  GMM: (converged)
```

### S3 methods

Standard model sumamry functions on the dbc output are provided.

Returning the bias-corrected coefficients:

``` r
coef(fit_endo)
#>          lag_y              D treat_eq:lag_y 
#>      0.1960087      0.4720423      0.2836029
```

Returning the variance-convariance matrix on the bias-corrected
coefficients:

``` r
vcov(fit_endo)
#>                        lag_y             D treat_eq:lag_y
#> lag_y           0.0005870888 -0.0001764293   0.0000619551
#> D              -0.0001764293  0.0009897852   0.0001792766
#> treat_eq:lag_y  0.0000619551  0.0001792766   0.0004988017
```

Standard confidence intervals:

``` r
confint(fit_endo, level = 0.90)
#>                       5%       95%
#> lag_y          0.1561540 0.2358634
#> D              0.4202938 0.5237908
#> treat_eq:lag_y 0.2468670 0.3203389
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
#> 1   Outcome lag_y 0.1960087 0.02422991  8.089533 5.989376e-16 0.1485189 0.2434985
#> 2   Outcome     D 0.4720423 0.03146085 15.004117 6.900337e-51 0.4103802 0.5337044
#> 3 Treatment lag_y 0.2836029 0.02233387 12.698333 6.040256e-37 0.2398293 0.3273765
```

`glance` shows other model information, like number of observations.

``` r
glance(fit_endo)
#>         N N_T nobs  phi converged
#> lag_y 500   3 1500 0.33      TRUE
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
| Outcome   | lag_y    | 0.215     | 0.196      |
|           |          | (0.021)   | (0.024)    |
|           | D        | 0.464     | 0.472      |
|           |          | (0.031)   | (0.031)    |
|           | X1_1     | 0.952     |            |
|           |          | (0.034)   |            |
| Treatment | lag_y    |           | 0.284      |
|           |          |           | (0.022)    |
|           | Num.Obs. | 1500      | 1500       |
|           | N        | 500       | 500        |
|           | N_T      | 3         | 3          |
|           | phi      | 0.215     | 0.33       |
