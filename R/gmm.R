#' Compute φ = ρ₁ + Σ_c τ_c ρ₂ μ_c
#'
#' Simple model (below Lemma 4.1): φ = ρ₁ + τ·ρ₂
#' Interaction model (below eq 104): φ = ρ₁ + Σ_c τ_c ρ₂ μ_c
#'
#' @keywords internal
compute_phi <- function(theta, parsed, mu_W, idx) {
    rho1 <- theta[idx$rho1]

    if (!parsed$has_treatment_eq || !parsed$has_lag_y_in_treatment) {
        return(rho1)
    }

    rho2 <- theta[idx$rho2]
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


#' Build an index mapping into the θ vector
#'
#' θ ordering: [ρ₁, τ_1,...,τ_{Nc}, β₁_1,..., ρ₂, β₂_1,...]
#'
#' @keywords internal
build_theta_index <- function(parsed) {
    pos <- 1

    rho1 <- pos
    pos <- pos + 1

    n_tau <- parsed$n_tau
    tau <- seq(pos, length.out = n_tau)
    pos <- pos + n_tau

    n_beta1 <- parsed$n_beta1
    beta1 <- if (n_beta1 > 0) seq(pos, length.out = n_beta1) else integer(0)
    pos <- pos + n_beta1

    rho2 <- integer(0)
    beta2 <- integer(0)

    if (parsed$has_treatment_eq) {
        if (parsed$has_lag_y_in_treatment) {
            rho2 <- pos
            pos <- pos + 1
        }
        n_beta2 <- parsed$n_beta2
        if (n_beta2 > 0) {
            beta2 <- seq(pos, length.out = n_beta2)
            pos <- pos + n_beta2
        }
    }

    list(
        rho1 = rho1,
        tau = tau,
        beta1 = beta1,
        rho2 = rho2,
        beta2 = beta2,
        n_params = pos - 1
    )
}


#' Kiviet-style bias kernel (for b_{ρ₁} and b_{τ})
#'
#' K(σ², φ, T) = -σ²/T² · [(T-1)/(1-φ) - (φ - φ^T)/(1-φ)²]
#'
#' @param sigma_sq Numeric vector: estimated variance (per-unit or scalar)
#' @param phi Scalar: φ(θ)
#' @param N_T Integer: number of time periods T
#' @return Numeric: bias term (same length as sigma_sq)
#' @keywords internal
kiviet_kernel <- function(sigma_sq, phi, N_T) {
    -sigma_sq / N_T^2 * ((N_T - 1) / (1 - phi) - (phi - phi^N_T) / (1 - phi)^2)
}


#' Kiviet-style bias kernel for b_{ρ₂}
#'
#' @keywords internal
kiviet_kernel_rho2 <- function(sigma_sq, phi, N_T) {
    # Same as kiviet_kernel
    kiviet_kernel(sigma_sq, phi, N_T)
}


#' GMM moment function for DBC
#'
#' Implements bias-corrected moment conditions from eq (32)-(33):
#'   m^{DBC}_{iT}(θ) = m_{iT}(θ) - b̂_{iT}(θ)
#'
#' Returns an N × n_moments matrix (one row per unit, time-averaged
#' per eq 34).
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
    rho1 <- theta[idx$rho1]
    tau <- theta[idx$tau]
    beta1 <- if (length(idx$beta1) > 0) theta[idx$beta1] else numeric(0)
    rho2 <- if (length(idx$rho2) > 0) theta[idx$rho2] else 0
    beta2 <- if (length(idx$beta2) > 0) theta[idx$beta2] else numeric(0)

    # Compute phi
    phi <- compute_phi(theta, parsed, mu_W, idx)

    # Outcome residuals ε̃(θ)
    # ε̃ = ỹ - X̃_out %*% β_out  where β_out = [ρ₁, τ, β₁]
    beta_out <- c(rho1, tau, beta1)
    # Reorder to match design matrix column order
    out_col_order <- c(dm_mats$rho1_cols, dm_mats$tau_cols, dm_mats$beta1_cols)
    eps_tilde <- y_dm - X_out[, out_col_order, drop = FALSE] %*% beta_out

    # Treatment residuals ũ(θ) (if treatment equation present)
    u_tilde <- NULL
    if (parsed$has_treatment_eq && !is.null(X_treat)) {
        beta_treat <- c(
            if (length(idx$rho2) > 0) rho2 else NULL,
            beta2
        )
        treat_col_order <- c(dm_mats$rho2_cols, dm_mats$beta2_cols)
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

    # Bias corrections per unit (eqs 29-31, 105)
    # b_{ρ₁,iT}(θ) = K(σ²_{ε,i}, φ, T)          (eq 29)
    b_rho1_i <- kiviet_kernel(sigma_eps_i, phi, N_T)

    # Build moment conditions
    n_moments <- idx$n_params
    n_obs <- length(eps_tilde)
    moments <- matrix(0, nrow = n_obs, ncol = n_moments)
    col <- 1

    # m^{DBC}_{ρ₁} = Ỹ_{t-1} · ε̃ - b_{ρ₁,i}   (eq 33, row 1)
    Y_lag_dm <- X_out[, dm_mats$rho1_cols]
    moments[, col] <- Y_lag_dm * eps_tilde - b_rho1_i
    col <- col + 1

    # τ moments (eq 33 row 2 / eq 105 row 2)
    tau_instruments <- X_out[, dm_mats$tau_cols, drop = FALSE]
    # Bias for τ terms:
    #   Simple: b_{τ} = ρ₂ · b_{ρ₁}                    (eq 30)
    #   Interaction: b_{τ_c} = ρ₂ · μ_c · K(σ²_ε, φ, T) (eq 105, row 2)
    mu_vec <- c(
        if (parsed$has_bare_treatment) 1 else NULL,
        mu_W
    )
    for (k in seq_along(dm_mats$tau_cols)) {
        b_tau_k <- rho2 * mu_vec[k] * kiviet_kernel(sigma_eps_i, phi, N_T)
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

    # Treatment equation moments (if present)
    if (parsed$has_treatment_eq && !is.null(u_tilde)) {
        if (parsed$has_lag_y_in_treatment) {
            # m^{DBC}_{ρ₂} = Ỹ_{t-1} · ũ - b_{ρ₂,i}  (eq 33, row 3)
            # b_{ρ₂,iT} = τ · K₂(σ²_{u,i}, φ, T)  (eq 31, note T² denominator)
            # Interaction model: b_{ρ₂} = (Σ τ_c μ_c) · K₂(σ²_{u,i}, φ, T) (eq 105 row 3)
            sum_tau_mu <- sum(tau * mu_vec)
            b_rho2_i <- sum_tau_mu * kiviet_kernel_rho2(sigma_u_i, phi, N_T)
            Y_lag_treat <- X_treat[, dm_mats$rho2_cols]
            moments[, col] <- Y_lag_treat * u_tilde - b_rho2_i
            col <- col + 1
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
