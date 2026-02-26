library(data.table)
library(gmm)
library(plm)
library(testthat)
library(tidyverse)

source("R/helpers.R")
source("R/gmm.R")
source("R/dgp.R")
source("R/dbc.R")

set.seed(12345)

# Generator data with lag Y in treatment formula
data <- DGP(N = 500, N_T = 4, rho1 = 0.2, rho2 = 0.3, tau = 0.5)

result <- dbc(
    outcome_fml = y ~ lag_y + D,
    treatment_fml = D ~ lag_y,
    lag_y = "lag_y",
    treatment = "D",
    panel_id = "panel_id",
    time_id = "time",
    data = data,
    reltol = 1e-12
)

# Test the convergence of the debiasing
# Test gmm convergence code is 0
expect_equal(result$gmm_fit$algoInfo$convergence, 0)
# Test that debiased estimates are close enough to true
expect_equal(unname(abs(result$coefficients[1] - 0.2) < 1e-2), TRUE)
expect_equal(unname(abs(result$coefficients[2] - 0.5) < 1e-2), TRUE)
expect_equal(unname(abs(result$coefficients[3] - 0.3) < 1e-2), TRUE)
# Did we just get three parameter estimates back?
expect_equal(length(result$coefficients), 3)

### Monte carlo test
set.seed(12345)

# Do 250 tests
N <- 250
test_results <- map(
    1:N,
    \(n) {
        data <- DGP(N = 500, N_T = 4, rho1 = 0.2, rho2 = 0.3, tau = 0.5)
        result <- dbc(
            outcome_fml = y ~ lag_y + D,
            treatment_fml = D ~ lag_y,
            lag_y = "lag_y",
            treatment = "D",
            panel_id = "panel_id",
            time_id = "time",
            data = data,
            reltol = 1e-12
        )
        return(result)
    }
)

# Get the biased and corrected results
coef_results <- map_df(
    test_results,
    \(x) {
        data.table(
            dbc_lag_y = x$coefficients[1],
            dbc_D = x$coefficients[2],
            dbc_treatmentlag_y = x$coefficients[3],
            biased_lag_y = x$biased_coefficients[1],
            biased_D = x$biased_coefficients[2],
            biased_treatmentlag_y = x$biased_coefficients[3]
        )
    }
)

# Calculate mse
coef_results %>%
    summarise(
        dbc_lag_y = mean((dbc_lag_y - 0.2)^2),
        dbc_D = mean((dbc_D - 0.5)^2),
        dbc_treatmentlag_y = mean((dbc_treatmentlag_y - 0.3)^2),
        biased_lag_y = mean((biased_lag_y - 0.2)^2),
        biased_D = mean((biased_D - 0.5)^2),
        biased_treatmentlag_y = mean((biased_treatmentlag_y - 0.3)^2)
    ) %>%
    pivot_longer(
        cols = everything(),
        names_to = c("estimator", "parameter"),
        names_sep = "_",
        values_to = "mse"
    ) %>%
    ggplot() +
    geom_col(
        aes(x = parameter, y = mse, fill = estimator),
        position = "dodge"
    ) +
    labs(x = "Coefficient", y = "MSE") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_discrete(name = NULL)
