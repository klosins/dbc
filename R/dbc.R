#' Dynamic Bias Correction (DBC) Estimator
#'
#' Implements the bias-corrected GMM estimator from Klosin (2024) for
#' dynamic panel models with fixed effects. Corrects for both Nickell bias
#' and dynamic bias.
#'
#' @param outcome_fml Formula for the outcome equation. No fixed effects (no |).
#'   Example: y ~ lag_y1 + D + X1 + D:W1
#' @param treatment_fml Optional formula for the treatment equation. If NULL,
#'   treatment is assumed exogenous (\eqn{\rho_2 = 0}{rho_2 = 0}). Example: D ~ lag_y1 + X2
#' @param lag_y Character vector of lagged outcome variable names in the outcome
#'   equation, ordered lag-1, lag-2, ... Currently only length-1 is supported.
#'   Example: "lag_y1", or c("lag_y1", "lag_y2") (latter errors until implemented).
#' @param lag_d Character vector of lagged treatment variable names in the
#'   treatment equation. Currently not supported (errors if non-empty).
#' @param treatment Character: name of the treatment variable.
#' @param panel_id Character: name of the panel unit identifier.
#' @param time_id Character: name of the time variable.
#' @param data A data.frame or data.table.
#' @param balance_type Character: passed to plm::make.pbalanced(). Default "shared.times".
#' @param reltol Double: Relative convergence tolerance of GMM objective function.
#' @param maxit Double: Maximum number of solver iterations.
#'
#' @return An object of class "dbc".
#' @examples
#' set.seed(1)
#'
#' # Exogenous treatment
#' dat_exog <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0, tau = 0.5)
#' fit_exog <- dbc(
#'     outcome_fml = y ~ lag_y + D,
#'     lag_y       = "lag_y",
#'     treatment   = "D",
#'     panel_id    = "panel_id",
#'     time_id     = "time",
#'     data        = dat_exog
#' )
#' print(fit_exog)
#'
#' # Endogenous treatment
#' dat_endo <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit_endo <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y         = "lag_y",
#'     treatment     = "D",
#'     panel_id      = "panel_id",
#'     time_id       = "time",
#'     data          = dat_endo
#' )
#' print(fit_endo)
#'
#' # Heterogeneous treatment effects (interactions)
#' dat_int <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2,
#'                tau = c(0.8, 1.2), n_W = 2)
#' fit_int <- dbc(
#'     outcome_fml   = y ~ lag_y + D:W_1 + D:W_2,
#'     treatment_fml = D ~ lag_y,
#'     lag_y         = "lag_y",
#'     treatment     = "D",
#'     panel_id      = "panel_id",
#'     time_id       = "time",
#'     data          = dat_int
#' )
#' print(fit_int)
#' @export
dbc <- function(
    outcome_fml,
    treatment_fml = NULL,
    lag_y,
    lag_d = character(0),
    treatment,
    panel_id,
    time_id,
    data,
    balance_type = "shared.times",
    reltol = 1e-12,
    maxit = 10000
) {
    # Convert to data.table
    if (nrow(data) > 0) {
        data <- data.table::as.data.table(data)
    } else {
        stop("Data has 0 rows.")
    }

    # Validate formulae, etc are valid inputs
    validate_inputs(
        outcome_fml,
        treatment_fml,
        lag_y,
        lag_d,
        treatment,
        panel_id,
        time_id,
        data
    )

    # Balance the panel
    data <- balance_panel(data, panel_id, time_id, balance_type)

    # Validate whether treatment is binary
    n_treat_vals <- data.table::uniqueN(data[[treatment]])
    if (n_treat_vals <= 2) {
        warning(
            "Treatment '",
            treatment,
            "' has only ",
            n_treat_vals,
            " unique value(s). The DBC bias formula assumes time-varying ",
            "treatment. Absorbing/binary treatment has a different (smaller) bias."
        )
    }

    # Find # of panel ids
    N <- data.table::uniqueN(data[[panel_id]])

    # Find # of time periods
    N_T <- data.table::uniqueN(data[[time_id]])

    # Parse the formulae
    parsed <- parse_formulas(
        outcome_fml,
        treatment_fml,
        lag_y,
        lag_d,
        treatment,
        data
    )

    # Demean data by panel id
    dm_mats <- build_demeaned_matrices(
        data,
        parsed,
        lag_y,
        lag_d,
        treatment,
        panel_id,
        time_id
    )

    # Solve OLS - returns biased values
    # Use these as debiasing starting values
    biased_ols <- get_biased_ols(dm_mats, parsed)

    # Solve for debiased results
    gmm_result <- run_dbc_gmm(
        dm_mats = dm_mats,
        parsed = parsed,
        N = N,
        N_T = N_T,
        start_vals = biased_ols,
        reltol = reltol,
        maxit = maxit
    )

    # Extract out and name the debiased results
    theta_hat <- gmm_result$coefficients
    names(theta_hat) <- names(biased_ols)

    # Check stationarity assumptions (equation 35)
    idx <- build_theta_index(parsed)
    phi_hat <- compute_phi(theta_hat, parsed, dm_mats$mu_W, idx)
    if (abs(phi_hat) >= 1) {
        warning(
            "Estimated |phi| = ",
            round(abs(phi_hat), 4),
            " >= 1. Stationarity assumption (Assumption 5, eq 35) violated. ",
            "Bias correction may be unreliable."
        )
    }

    # Return
    structure(
        list(
            coefficients = theta_hat,
            biased_coefficients = biased_ols,
            vcov = gmm_result$vcov,
            se = sqrt(diag(gmm_result$vcov)),
            phi = phi_hat,
            N = N,
            N_T = N_T,
            gmm_fit = gmm_result$fit,
            parsed = parsed,
            call = match.call()
        ),
        class = "dbc"
    )
}
