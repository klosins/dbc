#' Generate data from dynamic panel (Klosin eqs 88-89)
#'
#' If n_W = 0 (no interactions):
#'   Y_{i,t} = a_i + τ₀ D_{i,t} + β₁'X₁_i + ρ₁ Y_{i,t-1} + ε_{i,t}
#'
#' If n_W > 0 (with interactions, eq 88):
#'   Y_{i,t} = a_i + Σ_c τ_c (D_{i,t} · W^c_i) + β₁'X₁_i + ρ₁ Y_{i,t-1} + ε_{i,t}
#'
#' Treatment equation (eq 89):
#'   D_{i,t} = c_i + ρ₂ Y_{i,t-1} + β₂'X₂_i + u_{i,t}
#'
#' @param N Number of individuals
#' @param N_T Number of time periods
#' @param rho1 Coefficient on lagged outcome in outcome equation (ρ₁)
#' @param rho2 Coefficient on lagged outcome in treatment equation (ρ₂)
#' @param tau Treatment effect. Scalar if n_W=0, vector of length n_W if n_W>0. Default 1.
#' @param n_X1 Number of exogenous covariates in outcome equation (default 0)
#' @param n_X2 Number of exogenous covariates in treatment equation (default 0)
#' @param n_W Number of covariates interacting with treatment (default 0)
#' @param beta_X1 Coefficients on X1 covariates (default NULL → vector of 1s)
#' @param beta_X2 Coefficients on X2 covariates (default NULL → vector of 1s)
#' @return data.table with columns: panel_id, time, a, D, y, lag_y, plus covariate columns
#' @export
DGP <- function(
    N = 100,
    N_T = 10,
    rho1 = 0.8,
    rho2 = 0,
    tau = 1,
    n_X1 = 0,
    n_X2 = 0,
    n_W = 0,
    beta_X1 = NULL,
    beta_X2 = NULL
) {
    #  Validate and set defaults for coefficients
    if (n_W > 0) {
        if (length(tau) == 1) {
            tau <- rep(tau, n_W)
        }
        stopifnot(length(tau) == n_W)
    } else {
        stopifnot(length(tau) == 1)
    }
    if (is.null(beta_X1)) {
        beta_X1 <- rep(1, n_X1)
    }
    if (is.null(beta_X2)) {
        beta_X2 <- rep(1, n_X2)
    }
    stopifnot(length(beta_X1) == n_X1)
    stopifnot(length(beta_X2) == n_X2)

    #  Draw fixed effects and initial conditions
    a <- rnorm(N, 0, 5)
    D0 <- rnorm(N, a, 1)

    #  Draw time-varying covariates (iid N(0,1) per individual-time)
    n_obs <- N * N_T
    # X1: outcome equation covariates
    X1_mat <- NULL
    if (n_X1 > 0) {
        X1_mat <- matrix(rnorm(n_obs * n_X1), nrow = n_obs, ncol = n_X1)
    }
    # X2: treatment equation covariates
    X2_mat <- NULL
    if (n_X2 > 0) {
        X2_mat <- matrix(rnorm(n_obs * n_X2), nrow = n_obs, ncol = n_X2)
    }
    # W: interaction covariates
    W_mat <- NULL
    if (n_W > 0) {
        W_mat <- matrix(rnorm(n_obs * n_W), nrow = n_obs, ncol = n_W)
    }

    #  Build base data.table
    panel_id <- rep(1:N, each = N_T)
    time <- rep(1:N_T, N)

    DT <- data.table::data.table(
        panel_id = panel_id,
        time = time,
        a = rep(a, each = N_T)
    )

    # Add covariate columns (time-varying, already in panel order)
    if (n_X1 > 0) {
        for (k in seq_len(n_X1)) {
            col_name <- paste0("X1_", k)
            DT[, (col_name) := X1_mat[, k]]
        }
    }
    if (n_X2 > 0) {
        for (k in seq_len(n_X2)) {
            col_name <- paste0("X2_", k)
            DT[, (col_name) := X2_mat[, k]]
        }
    }
    if (n_W > 0) {
        for (k in seq_len(n_W)) {
            col_name <- paste0("W_", k)
            DT[, (col_name) := W_mat[, k]]
        }
    }

    #  Helper to get row indices for a given time period
    # Panel is ordered: (id=1,t=1), (id=1,t=2), ..., (id=1,t=N_T), (id=2,t=1), ...
    rows_at <- function(t) seq(t, n_obs, by = N_T)

    #  Simulate panel
    # t=1: use initial D0
    r1 <- rows_at(1)

    # X2 contribution for treatment eq at t=1 (not used for D0, but compute for consistency)
    X2_contrib_t <- rep(0, N)
    if (n_X2 > 0) {
        X2_contrib_t <- as.numeric(X2_mat[r1, , drop = FALSE] %*% beta_X2)
    }

    # X1 contribution for outcome eq at t=1
    X1_contrib_t <- rep(0, N)
    if (n_X1 > 0) {
        X1_contrib_t <- as.numeric(X1_mat[r1, , drop = FALSE] %*% beta_X1)
    }

    DT[time == 1, D := rep(D0, each = 1)]

    if (n_W > 0) {
        W_t <- W_mat[r1, , drop = FALSE]
        DT[
            time == 1,
            y := a +
                {
                    interaction_sum <- 0
                    for (cc in seq_len(n_W)) {
                        interaction_sum <- interaction_sum +
                            tau[cc] * D * W_t[, cc]
                    }
                    interaction_sum
                } +
                X1_contrib_t +
                rnorm(.N, 0, 1)
        ]
    } else {
        DT[time == 1, y := a + tau * D + X1_contrib_t + rnorm(.N, 0, 1)]
    }

    for (year in 2:N_T) {
        DT[, lag_y := data.table::shift(y), by = panel_id]

        rt <- rows_at(year)

        X2_contrib_t <- rep(0, N)
        if (n_X2 > 0) {
            X2_contrib_t <- as.numeric(X2_mat[rt, , drop = FALSE] %*% beta_X2)
        }

        X1_contrib_t <- rep(0, N)
        if (n_X1 > 0) {
            X1_contrib_t <- as.numeric(X1_mat[rt, , drop = FALSE] %*% beta_X1)
        }

        # D_{i,t} = c_i + ρ₂ Y_{i,t-1} + β₂'X₂_{i,t} + u_{i,t}
        DT[time == year, D := a + rho2 * lag_y + X2_contrib_t + rnorm(.N, 0, 1)]

        if (n_W > 0) {
            W_t <- W_mat[rt, , drop = FALSE]
            DT[
                time == year,
                y := a +
                    {
                        interaction_sum <- 0
                        for (cc in seq_len(n_W)) {
                            interaction_sum <- interaction_sum +
                                tau[cc] * D * W_t[, cc]
                        }
                        interaction_sum
                    } +
                    X1_contrib_t +
                    rho1 * lag_y +
                    rnorm(.N, 0, 1)
            ]
        } else {
            DT[
                time == year,
                y := a + tau * D + X1_contrib_t + rho1 * lag_y + rnorm(.N, 0, 1)
            ]
        }
    }

    # Recompute lag_y for final dataset
    DT[, lag_y := data.table::shift(y), by = panel_id]

    # Drop first period (no lag available)
    DT <- DT[time > 1]

    DT
}
