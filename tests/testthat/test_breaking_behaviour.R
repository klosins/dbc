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

### Test giving variable names that don't exist
# Outcome formula
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D:W_3 + D:W_2 + X1_2,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data,
        reltol = 1e-12
    )
)

# Treatment formula
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
        treatment_fml = D ~ lag_y + X2_1 + AA,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data,
        reltol = 1e-12
    )
)

# Lag y
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
        treatment_fml = D ~ lag_y + X2_1 + AA,
        lag_y = "lag_y_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data,
        reltol = 1e-12
    )
)

# Treatment
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "DD",
        panel_id = "panel_id",
        time_id = "time",
        data = data,
        reltol = 1e-12
    )
)

# Panel id
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id_id",
        time_id = "time",
        data = data,
        reltol = 1e-12
    )
)

# Time
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time_time",
        data = data,
        reltol = 1e-12
    )
)

# Treatment not in outcome formula
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + X2_1:W_1 + X2_1:W_2 + X1_1,
        treatment_fml = D ~ lag_y,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data,
        reltol = 1e-12
    )
)

# Outcome variable on RHS
expect_error(
    dbc(
        outcome_fml = y ~ y + lag_y + D:W_1 + D:W_2 + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data,
        reltol = 1e-12
    )
)


### Data problems
# Given a binary treatment
expect_warning(
    dbc(
        outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data %>%
            mutate(D = sample(c(0, 1), size = nrow(data), replace = TRUE)),
        reltol = 1e-12
    )
)

# When phi is >= 1, stationarity assumptions are violated
# Generate data with high phi
data_high_phi <- DGP(
    N = 500,
    N_T = 4,
    rho1 = 0.9,
    rho2 = 0.9,
    tau = 0.8,
    n_X1 = 1,
    n_X2 = 1
)

# Get the warning
expect_warning(
    dbc(
        outcome_fml = y ~ lag_y + D + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data_high_phi,
        reltol = 1e-12
    )
)

# Data with no rows
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = data[NULL, ],
        balance_type = "shared.times",
        reltol = 1e-12
    )
)


# Unbalanced panel, do we get the warning?
data <- as.matrix(data[
    sample(1:nrow(data), floor(nrow(data) * 0.7)),
])
# Get the warning if we get positive data rows
expect_warning(
    dbc(
        outcome_fml = y ~ lag_y + D + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = unbalanced_data,
        balance_type = "shared.individuals",
        reltol = 1e-12
    )
)

# Get the stop error if we remove all rows
expect_error(
    dbc(
        outcome_fml = y ~ lag_y + D + X1_1,
        treatment_fml = D ~ lag_y + X2_1,
        lag_y = "lag_y",
        treatment = "D",
        panel_id = "panel_id",
        time_id = "time",
        data = unbalanced_data,
        balance_type = "shared.times",
        reltol = 1e-12
    )
)
