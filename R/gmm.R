#' Compute \eqn{\phi = \rho_1 + \sum_c \tau_c \rho_2 \mu_c}{phi = rho_1 + sum_c(tau_c*rho_2*mu_c)}
#'
#' Simple model (below Lemma 4.1): \eqn{\phi = \rho_1 + \tau \cdot \rho_2}{phi = rho_1 + tau*rho_2}
#'
#' Interaction model (below eq 104): \eqn{\phi = \rho_1 + \sum_c \tau_c \rho_2 \mu_c}{phi = rho_1 + sum_c(tau_c*rho_2*mu_c)}
#'
#' Currently only the first Y-lag coefficient (\eqn{\rho_1}{rho_1}) enters \eqn{\phi}{phi}.
#' When multi-lag bias corrections are derived, \eqn{\phi}{phi} will generalise.
#'
#' @keywords internal
compute_phi <- function(theta, parsed, mu_W, idx) {
    # φ is defined via the first Y-lag only (current supported case: n_lags_Y == 1)
    rho1 <- theta[idx$rho[1]]

    if (!parsed$has_treatment_eq || !parsed$has_lag_y_in_treatment) {
        return(rho1)
    }

    rho2 <- theta[idx$rho2[1]]
    tau <- theta[idx$tau]

    if (length(parsed$W_vars) > 0) {
        # Interaction model: φ = ρ₁ + Σ_c τ_c ρ₂ μ_c
        # tau vector aligns with: [bare D (if present), D:W1, D:W2, ...]
        mu_vec <- c(
            if (parsed$has_bare_treatment) 1 else NULL,
            mu_W
        )
        phi <- rho1 + rho2 * sum(tau * mu_vec)
    } else {
        # Simple model: φ = ρ₁ + τ·ρ₂
        phi <- rho1 + tau * rho2
    }
    phi
}


#' Build an index mapping into the \eqn{\theta}{theta} vector
#'
#' \eqn{\theta}{theta} ordering: \eqn{(\rho_1, \ldots, \rho_{n_L},\ \tau_1, \ldots, \tau_{N_c},\ \beta_{1,1}, \ldots,\ \rho_{2,1}, \ldots,\ \mathrm{lag\_d}_1, \ldots,\ \beta_{2,1}, \ldots)}{(rho_1,...,rho_nL, tau_1,...,tau_Nc, beta1_1,..., rho2_1,..., lag_d_1,..., beta2_1,...)}
#'
#' rho: integer vector of length n_lags_Y (indices of Y-lag coefficients)
#' rho2: integer vector of length n_rho2 (Y-lag coefs in treatment equation)
#' lag_d: integer vector of length n_lag_d (D-lag coefs in treatment equation, placeholder)
#'
#' @keywords internal
build_theta_index <- function(parsed) {
    pos <- 1

    # Y-lag coefficients in outcome equation: one per lag
    n_lags_Y <- parsed$n_lags_Y
    rho <- seq(pos, length.out = n_lags_Y)
    pos <- pos + n_lags_Y

    n_tau <- parsed$n_tau
    tau <- seq(pos, length.out = n_tau)
    pos <- pos + n_tau

    n_beta1 <- parsed$n_beta1
    beta1 <- if (n_beta1 > 0) seq(pos, length.out = n_beta1) else integer(0)
    pos <- pos + n_beta1

    rho2 <- integer(0)
    lag_d <- integer(0)
    beta2 <- integer(0)

    if (parsed$has_treatment_eq) {
        # Y-lag coefficients in treatment equation
        n_rho2 <- parsed$n_rho2
        if (n_rho2 > 0) {
            rho2 <- seq(pos, length.out = n_rho2)
            pos <- pos + n_rho2
        }

        # D-lag coefficients in treatment equation (placeholder for future bias terms)
        n_lag_d <- parsed$n_lag_d
        if (n_lag_d > 0) {
            lag_d <- seq(pos, length.out = n_lag_d)
            pos <- pos + n_lag_d
        }

        n_beta2 <- parsed$n_beta2
        if (n_beta2 > 0) {
            beta2 <- seq(pos, length.out = n_beta2)
            pos <- pos + n_beta2
        }
    }

    list(
        rho = rho, # length = n_lags_Y
        tau = tau,
        beta1 = beta1,
        rho2 = rho2, # length = n_rho2
        lag_d = lag_d, # length = n_lag_d (placeholder)
        beta2 = beta2,
        n_params = pos - 1
    )
}


#' Bias term core
#'
#' \deqn{-\frac{\sigma^2}{T^2} \left(\frac{T-1}{1-\phi} - \frac{\phi - \phi^T}{(1-\phi)^2}\right)}{-sigma^2/T^2 * ((T-1)/(1-phi) - (phi - phi^T)/(1-phi)^2)}
#'
#' @param sigma_sq Numeric vector: estimated variance (per-unit or scalar)
#' @param phi Scalar: \eqn{\phi(\theta)}{phi(theta)}
#' @param N_T Integer: number of time periods T
#' @return Numeric: bias term (same length as sigma_sq)
#' @keywords internal
bias_term_core <- function(sigma_sq, phi, N_T) {
    -sigma_sq / N_T^2 * ((N_T - 1) / (1 - phi) - (phi - phi^N_T) / (1 - phi)^2)
}


#' GMM moment function for DBC
#'
#' Implements bias-corrected moment conditions from eq (32)-(33):
#'
#' \deqn{m^{DBC}_{iT}(\theta) = m_{iT}(\theta) - \hat{b}_{iT}(\theta)}{m_DBC_iT(theta) = m_iT(theta) - b_hat_iT(theta)}
#'
#' Returns an \eqn{N \times n_{\mathrm{moments}}}{N x n_moments} matrix (one row per unit,
#' time-averaged per eq 34).
#'
#' @param theta Parameter vector
#' @param x List containing demeaned matrices and metadata
#' @return N × n_moments matrix
#' @keywords internal
gmm_moment_fn <- function(theta, x) {
    dm_mats <- x$dm_mats
    parsed <- x$parsed
    idx <- x$idx
    mu_W <- dm_mats$mu_W
    N <- x$N
    N_T <- x$N_T
    panel_vec <- dm_mats$panel_vec

    X_out <- dm_mats$X_out_dm
    X_treat <- dm_mats$X_treat_dm
    y_dm <- dm_mats$y_dm
    D_dm <- dm_mats$D_dm

    # Extract parameters
    # rho: vector of Y-lag coefficients [ρ₁, ρ₂, ...] in outcome equation
    rho <- theta[idx$rho]
    tau <- theta[idx$tau]
    beta1 <- if (length(idx$beta1) > 0) theta[idx$beta1] else numeric(0)
    # rho2: Y-lag coefficients in treatment equation [ρ₂₁, ρ₂₂, ...]
    rho2 <- if (length(idx$rho2) > 0) theta[idx$rho2] else numeric(0)
    # lag_d: D-lag coefficients in treatment equation (placeholder)
    lag_d_coef <- if (length(idx$lag_d) > 0) theta[idx$lag_d] else numeric(0)
    # Other explanatory variables in treatment equation
    beta2 <- if (length(idx$beta2) > 0) theta[idx$beta2] else numeric(0)

    # Compute phi
    phi <- compute_phi(theta, parsed, mu_W, idx)

    # Outcome residuals ε̃(θ)
    # ε̃ = ỹ - X̃_out %*% β_out  where β_out = [ρ₁,..., τ, β₁]
    beta_out <- c(rho, tau, beta1)
    # Reorder to match design matrix column order
    out_col_order <- c(
        unlist(dm_mats$rho_cols),
        dm_mats$tau_cols,
        dm_mats$beta1_cols
    )
    eps_tilde <- y_dm - X_out[, out_col_order, drop = FALSE] %*% beta_out

    # Treatment residuals ũ(θ) (if treatment equation present)
    u_tilde <- NULL
    if (parsed$has_treatment_eq && !is.null(X_treat)) {
        # β_treat = [ρ₂₁,..., lag_d coefs, β₂]
        beta_treat <- c(rho2, beta2)
        treat_col_order <- c(
            unlist(dm_mats$rho2_cols),
            unlist(dm_mats$lag_d_cols),
            dm_mats$beta2_cols
        )
        u_tilde <- D_dm -
            X_treat[, treat_col_order, drop = FALSE] %*% beta_treat
    }

    eps_tilde <- as.numeric(eps_tilde)
    if (!is.null(u_tilde)) {
        u_tilde <- as.numeric(u_tilde)
    }

    # Unit-level variances (eqs 27-28)
    # σ²_{ε,iT}(θ) = 1/(T-1) Σ_t ε̃²_{i,t}     (eq 27)
    sigma_eps_i <- ave(eps_tilde^2, panel_vec, FUN = function(x) {
        sum(x) / (N_T - 1)
    })

    sigma_u_i <- NULL
    if (!is.null(u_tilde)) {
        # σ²_{u,iT}(θ) = 1/(T-1) Σ_t ũ²_{i,t}  (eq 28)
        sigma_u_i <- ave(u_tilde^2, panel_vec, FUN = function(x) {
            sum(x) / (N_T - 1)
        })
    }

    # Build moment conditions
    n_moments <- idx$n_params
    n_obs <- length(eps_tilde)
    moments <- matrix(0, nrow = n_obs, ncol = n_moments)
    col <- 1

    # --- Outcome equation moments ---

    # m^{DBC}_{ρ₁} = Ỹ_{t-1} · ε̃ - b_{ρ₁,i}   (eq 33, row 1)
    # Currently only ρ₁ (first Y-lag) has a bias correction derived.
    Y_lag1_dm <- X_out[, dm_mats$rho_cols[[1]]]
    b_rho1_i <- bias_term_core(sigma_eps_i, phi, N_T)
    moments[, col] <- Y_lag1_dm * eps_tilde - b_rho1_i
    col <- col + 1

    # Higher Y-lag moments (no bias correction derived yet — plain OLS moment)
    if (parsed$n_lags_Y > 1) {
        for (k in 2:parsed$n_lags_Y) {
            # Not yet supported
            # Relevant data objects are:
            #   # New columnn in the moment matrix: moments[, col]
            #   # The Vector of each lag Y for k > 2: X_out[, dm_mats$rho_cols[[k]]]
            #   # Update moment columns: col <- col + 1
        }
    }

    # τ moments (eq 33 row 2 / eq 105 row 2)
    tau_instruments <- X_out[, dm_mats$tau_cols, drop = FALSE]
    # Bias for τ terms:
    #   Simple: b_{τ} = ρ₂ · b_{ρ₁}                    (eq 30)
    #   Interaction: b_{τ_c} = ρ₂ · μ_c · K(σ²_ε, φ, T) (eq 105, row 2)
    mu_vec <- c(
        if (parsed$has_bare_treatment) 1 else NULL,
        mu_W
    )
    rho2_scalar <- if (length(rho2) > 0) rho2[1] else 0
    for (k in seq_along(dm_mats$tau_cols)) {
        b_tau_k <- rho2_scalar *
            mu_vec[k] *
            bias_term_core(sigma_eps_i, phi, N_T)
        moments[, col] <- tau_instruments[, k] * eps_tilde - b_tau_k
        col <- col + 1
    }

    # β₁ moments: no bias correction (eq 105, row 4: b_{β₁} = 0)
    if (length(dm_mats$beta1_cols) > 0) {
        X1_dm <- X_out[, dm_mats$beta1_cols, drop = FALSE]
        for (k in seq_len(ncol(X1_dm))) {
            moments[, col] <- X1_dm[, k] * eps_tilde
            col <- col + 1
        }
    }

    # --- Treatment equation moments (if present) ---
    if (parsed$has_treatment_eq && !is.null(u_tilde)) {
        if (parsed$has_lag_y_in_treatment) {
            # m^{DBC}_{ρ₂} = Ỹ_{t-1} · ũ - b_{ρ₂,i}  (eq 33, row 3)
            # b_{ρ₂,iT} = τ · K₂(σ²_{u,i}, φ, T)  (eq 31, note T² denominator)
            # Interaction model: b_{ρ₂} = (Σ τ_c μ_c) · K₂(σ²_{u,i}, φ, T) (eq 105 row 3)
            sum_tau_mu <- sum(tau * mu_vec)
            b_rho2_i <- sum_tau_mu * bias_term_core(sigma_u_i, phi, N_T)
            Y_lag_treat <- X_treat[, dm_mats$rho2_cols[[1]]]
            moments[, col] <- Y_lag_treat * u_tilde - b_rho2_i
            col <- col + 1
        }

        # D-lag moments in treatment equation (placeholder)
        if (length(dm_mats$lag_d_cols) > 0) {
            for (lag_d_col_idx in unlist(dm_mats$lag_d_cols)) {
                # Not yet supported
                # Relevant data objects are:
                #   # New columnn in the moment matrix: moments[, col]
                #   # The Vector of each lag D:  X_treat[, lag_d_col_idx]
                #   # Update moment columns: col <- col + 1
            }
        }

        # β₂ moments: no bias correction (eq 105, row 5: b_{β₂} = 0)
        if (length(dm_mats$beta2_cols) > 0) {
            X2_dm <- X_treat[, dm_mats$beta2_cols, drop = FALSE]
            for (k in seq_len(ncol(X2_dm))) {
                moments[, col] <- X2_dm[, k] * u_tilde
                col <- col + 1
            }
        }
    }

    # Average moments within each unit over T (eq 34)
    # GMM objective (eq 34) is over m^{DBC}_{iT} = (1/T) Σ_t moments_{i,t}
    # gmm::gmm then computes (1/N) Σ_i internally.
    moments_dt <- data.table::as.data.table(moments)
    moments_dt[, .panel_grp := panel_vec]
    col_nms <- paste0("V", seq_len(n_moments))
    moments_agg <- moments_dt[,
        lapply(.SD, mean),
        by = .panel_grp,
        .SDcols = col_nms
    ]
    moments_agg[, .panel_grp := NULL]
    as.matrix(moments_agg)
}


#' Run the DBC GMM estimation
#' @keywords internal
run_dbc_gmm <- function(
    dm_mats,
    parsed,
    N,
    N_T,
    start_vals,
    reltol = 1e-12,
    maxit = 10000
) {
    idx <- build_theta_index(parsed)

    x_data <- list(
        dm_mats = dm_mats,
        parsed = parsed,
        idx = idx,
        N = N,
        N_T = N_T
    )

    # GMM estimation
    fit <- gmm::gmm(
        g = gmm_moment_fn,
        x = x_data,
        t0 = start_vals,
        type = "twoStep",
        wmatrix = "ident",
        vcov = "iid",
        method = "BFGS",
        control = list(reltol = reltol, maxit = maxit)
    )

    list(
        coefficients = coef(fit),
        vcov = vcov(fit),
        fit = fit
    )
}
