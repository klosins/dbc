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

# Monte Carlo simulation of interaction model

# Do 250 tests
N <- 250
test_results <- map(
    1:N,
    \(n) {
        # Generate data with lag Y in treatment formula
        data <- DGP(
            N = 500,
            N_T = 4,
            rho1 = 0.2,
            rho2 = 0.3,
            tau = c(0.8, 1.2),
            n_X1 = 1,
            n_X2 = 1,
            n_W = 2
        )
        result <- dbc(
            outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
            treatment_fml = D ~ lag_y + X2_1,
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
            dbc_DxW1 = x$coefficients[2],
            dbc_DxW2 = x$coefficients[3],
            dbc_treatmentlag_y = x$coefficients[5],
            biased_lag_y = x$biased_coefficients[1],
            biased_DxW1 = x$biased_coefficients[2],
            biased_DxW2 = x$biased_coefficients[3],
            biased_treatmentlag_y = x$biased_coefficients[5]
        )
    }
)

# Chart of MSE
g <-
    coef_results %>%
    summarise(
        `dbc.Lag Y` = mean((dbc_lag_y - 0.2)^2),
        dbc.DxW1 = mean((dbc_DxW1 - 0.8)^2),
        dbc.DxW2 = mean((dbc_DxW2 - 1.2)^2),
        `dbc.Treatment equation lag Y` = mean((dbc_treatmentlag_y - 0.3)^2),
        `Biased.Lag Y` = mean((biased_lag_y - 0.2)^2),
        Biased.DxW1 = mean((biased_DxW1 - 0.8)^2),
        Biased.DxW2 = mean((biased_DxW2 - 1.2)^2),
        `Biased.Treatment equation lag Y` = mean(
            (biased_treatmentlag_y - 0.3)^2
        )
    ) %>%
    pivot_longer(
        cols = everything(),
        names_to = c("estimator", "parameter"),
        names_sep = "\\.",
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
ggsave(
    "readme_code/interaction_model_mse.png",
    plot = g,
    width = 11,
    height = 5,
    units = "in"
)
