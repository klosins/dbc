#' Validate inputs to dbc()
#' @keywords internal
validate_inputs <- function(
    outcome_fml,
    treatment_fml,
    lag_y,
    lag_d,
    treatment,
    panel_id,
    time_id,
    data
) {
    stopifnot(
        "outcome_fml must be a formula" = inherits(outcome_fml, "formula"),
        "lag_y must be a character vector" = is.character(lag_y) &&
            length(lag_y) >= 1,
        "lag_d must be a character vector" = is.character(lag_d),
        "treatment must be a character string" = is.character(treatment) &&
            length(treatment) == 1,
        "panel_id must be a character string" = is.character(panel_id) &&
            length(panel_id) == 1,
        "time_id must be a character string" = is.character(time_id) &&
            length(time_id) == 1
    )

    # Multiple Y lags or any D lags not yet supported
    if (length(lag_y) > 1) {
        stop(
            "Multiple lags of Y (lag_y of length > 1) are not yet supported. ",
            "Bias correction for higher-order Y lags is not implemented."
        )
    }
    if (length(lag_d) > 0) {
        stop(
            "Lags of D in the treatment equation (lag_d) are not yet supported. ",
            "Bias correction for lagged treatment terms is not implemented."
        )
    }

    if (!is.null(treatment_fml)) {
        stopifnot(
            "treatment_fml must be a formula" = inherits(
                treatment_fml,
                "formula"
            )
        )
        treat_lhs <- all.vars(treatment_fml[[2]])
        if (!identical(treat_lhs, treatment)) {
            stop(
                "LHS of treatment_fml must be '",
                treatment,
                "', got '",
                paste(treat_lhs, collapse = ", "),
                "'"
            )
        }
    }

    # Check no FE syntax in formulas
    if (grepl("\\|", deparse(outcome_fml))) {
        stop(
            "outcome_fml should not contain '|'. ",
            "Demeaning is handled internally by panel_id. ",
            "Additional FEs should enter as covariates."
        )
    }
    if (!is.null(treatment_fml) && grepl("\\|", deparse(treatment_fml))) {
        stop(
            "treatment_fml should not contain '|'. ",
            "Demeaning is handled internally by panel_id. ",
            "Additional FEs should enter as covariates."
        )
    }

    # Check variables exist in data
    all_vars <- unique(c(all.vars(outcome_fml), lag_d, panel_id, time_id))
    if (!is.null(treatment_fml)) {
        all_vars <- unique(c(all_vars, all.vars(treatment_fml)))
    }
    missing_vars <- setdiff(all_vars, names(data))
    if (length(missing_vars) > 0) {
        stop(
            "Variables not found in data: ",
            paste(missing_vars, collapse = ", ")
        )
    }

    # Check all lag_y appear in outcome RHS
    outcome_rhs_vars <- all.vars(outcome_fml[[3]])
    missing_lag_y <- setdiff(lag_y, outcome_rhs_vars)
    if (length(missing_lag_y) > 0) {
        stop(
            "lag_y variable(s) not found on RHS of outcome_fml: ",
            paste(missing_lag_y, collapse = ", ")
        )
    }

    # Check treatment appears in outcome RHS (as main effect or in interaction)
    if (!treatment %in% outcome_rhs_vars) {
        stop(
            "treatment ('",
            treatment,
            "') must appear on the RHS of outcome_fml."
        )
    }
}


#' Balance panel using plm
#' @keywords internal
balance_panel <- function(data, panel_id, time_id, balance_type) {
    pdata <- plm::pdata.frame(data, index = c(panel_id, time_id))
    if (plm::is.pbalanced(pdata)) {
        return(data)
    }
    message(
        "Panel is unbalanced. Balancing with balance.type = '",
        balance_type,
        "'."
    )
    balanced <- plm::make.pbalanced(pdata, balance.type = balance_type)
    if (nrow(balanced) == 0) {
        stop(
            "Balancing with '",
            balance_type,
            "' removes all data rows."
        )
    }
    dt <- data.table::as.data.table(balanced)
    if (!panel_id %in% names(dt)) {
        idx <- attr(balanced, "index")
        dt[[panel_id]] <- idx[[1]]
        dt[[time_id]] <- idx[[2]]
    }

    dt
}


#' Parse formulas using terms() and the factors attribute
#'
#' Classifies each term in the outcome equation as:
#'   - rho_terms: lag_y main effects (one per lag, ordered lag-1, lag-2, ...)
#'   - tau: any term involving treatment (bare D, D:W interactions)
#'   - beta1: everything else (exogenous controls, their interactions, etc.)
#'
#' For the treatment equation:
#'   - rho2_terms: lag_y terms (one per Y-lag present in treatment equation)
#'   - lag_d_terms: lag_d terms (D lags in treatment equation, placeholder)
#'   - beta2: everything else
#'
#' Uses attr(terms(), "factors") — a matrix where rows = variables,
#' cols = terms, entries > 0 indicate involvement. No regex needed.
#'
#' @param lag_y Character vector of lagged outcome variable names,
#'   ordered lag-1, lag-2, ... (currently only length-1 is supported).
#' @param lag_d Character vector of lagged treatment variable names
#'   in the treatment equation (currently only length-0 is supported).
#' @keywords internal
parse_formulas <- function(
    outcome_fml,
    treatment_fml,
    lag_y,
    lag_d,
    treatment,
    data
) {
    y_var <- all.vars(outcome_fml[[2]])
    if (length(y_var) != 1) {
        stop("Outcome formula must have a single LHS variable.")
    }

    # Extract out the outcome equation
    outcome_terms <- terms(outcome_fml, data = data)
    outcome_labels <- attr(outcome_terms, "term.labels")
    outcome_factors <- attr(outcome_terms, "factors")

    # For each lag_y variable, find which outcome terms involve it
    # involves_any_lag_y[j] = TRUE if term j involves any Y lag
    involves_lag_y_by_var <- lapply(lag_y, function(ly) {
        if (ly %in% rownames(outcome_factors)) {
            outcome_factors[ly, ] > 0
        } else {
            rep(FALSE, length(outcome_labels))
        }
    })
    involves_any_lag_y <- Reduce(`|`, involves_lag_y_by_var)

    involves_treatment <- if (treatment %in% rownames(outcome_factors)) {
        outcome_factors[treatment, ] > 0
    } else {
        rep(FALSE, length(outcome_labels))
    }

    # Assign bias correcting terms to each of the coefficient types
    # rho_terms: one entry per lag_y variable (the main-effect term for that lag)
    is_tau <- involves_treatment & !involves_any_lag_y
    is_beta1 <- !involves_any_lag_y & !involves_treatment

    if (any(involves_any_lag_y & involves_treatment)) {
        stop("Interactions between lag_y and treatment are not supported.")
    }

    # Collect the term label for each Y-lag (must be a bare main effect)
    rho_terms <- character(length(lag_y))
    for (k in seq_along(lag_y)) {
        is_rho_k <- involves_lag_y_by_var[[k]] & !involves_treatment
        if (sum(is_rho_k) != 1) {
            stop(
                "Expected exactly one term for lag_y '",
                lag_y[k],
                "' in outcome equation, found ",
                sum(is_rho_k),
                if (sum(is_rho_k) > 0) {
                    paste0(
                        ": ",
                        paste(outcome_labels[is_rho_k], collapse = ", ")
                    )
                } else {
                    ""
                }
            )
        }
        rho_terms[k] <- outcome_labels[is_rho_k]
    }

    tau_terms <- outcome_labels[is_tau]
    beta1_terms <- outcome_labels[is_beta1]

    # Identify W variables from interaction terms
    W_vars <- character(0)
    has_bare_treatment <- FALSE
    for (tt in tau_terms) {
        term_vars <- names(which(outcome_factors[, tt] > 0))
        other_vars <- setdiff(term_vars, treatment)
        if (length(other_vars) == 0) {
            has_bare_treatment <- TRUE
        } else {
            W_vars <- c(W_vars, other_vars)
        }
    }

    # Do the same for the treatment equation
    has_treatment_eq <- !is.null(treatment_fml)
    rho2_terms <- character(0) # Y-lag terms in treatment equation
    lag_d_terms <- character(0) # D-lag terms in treatment equation (placeholder)
    beta2_terms <- character(0)

    if (has_treatment_eq) {
        treat_term_obj <- terms(treatment_fml, data = data)
        treat_labels <- attr(treat_term_obj, "term.labels")
        treat_factors <- attr(treat_term_obj, "factors")

        # Y-lag terms in treatment equation (ordered as lag_y)
        rho2_terms <- intersect(lag_y, treat_labels)

        # D-lag terms in treatment equation (ordered as lag_d)
        lag_d_terms <- intersect(lag_d, treat_labels)

        beta2_terms <- setdiff(treat_labels, c(rho2_terms, lag_d_terms))
    }

    has_lag_y_in_treatment <- length(rho2_terms) > 0

    list(
        y_var = y_var,
        outcome_fml = outcome_fml,
        treatment_fml = treatment_fml,
        # Outcome equation
        rho_terms = rho_terms, # length = n_lags_Y
        tau_terms = tau_terms,
        beta1_terms = beta1_terms,
        W_vars = W_vars,
        has_bare_treatment = has_bare_treatment,
        # Treatment equation
        has_treatment_eq = has_treatment_eq,
        has_lag_y_in_treatment = has_lag_y_in_treatment,
        rho2_terms = rho2_terms, # Y-lag terms in treatment eq
        lag_d_terms = lag_d_terms, # D-lag terms in treatment eq (placeholder)
        beta2_terms = beta2_terms,
        # Counts
        n_lags_Y = length(lag_y),
        n_lags_D = length(lag_d),
        n_tau = length(tau_terms),
        n_beta1 = length(beta1_terms),
        n_rho2 = length(rho2_terms),
        n_lag_d = length(lag_d_terms),
        n_beta2 = length(beta2_terms)
    )
}


#' Build design matrices via model.matrix() and demean by panel
#'
#' model.matrix() handles all interaction expansion, factor coding, etc.
#' We then within-transform (demean) by panel_id.
#'
#' @return List with demeaned matrices, vectors, column mappings, and μ_c.
#' @keywords internal
build_demeaned_matrices <- function(
    data,
    parsed,
    lag_y,
    lag_d,
    treatment,
    panel_id,
    time_id
) {
    # Outcome design matrix
    X_out_raw <- model.matrix(parsed$outcome_fml, data = data)
    X_out <- X_out_raw
    # Model matrix gives an intercept by default, so remove that
    if ("(Intercept)" %in% colnames(X_out)) {
        X_out <- X_out[, colnames(X_out) != "(Intercept)", drop = FALSE]
    }

    y_raw <- data[[parsed$y_var]]
    D_raw <- data[[treatment]]

    # Treatment design matrix
    X_treat <- NULL
    X_treat_raw <- NULL
    if (parsed$has_treatment_eq) {
        X_treat_raw <- model.matrix(parsed$treatment_fml, data = data)
        X_treat <- X_treat_raw
        # Again, drop the intercept
        if ("(Intercept)" %in% colnames(X_treat)) {
            X_treat <- X_treat[,
                colnames(X_treat) != "(Intercept)",
                drop = FALSE
            ]
        }
    }

    # For interaction models where we have DxW in the outcome equation
    # We need mu for debiasing
    # Mu is just the mean of W (equation 105)
    if (length(parsed$W_vars) > 0) {
        mu_W <- vapply(
            parsed$W_vars,
            function(w) mean(data[[w]], na.rm = TRUE),
            numeric(1)
        )
    } else {
        mu_W <- numeric(0)
    }

    # Demean by panel_id
    panel_vec <- data[[panel_id]]

    demean_vec <- function(x) {
        x - ave(x, panel_vec, FUN = mean)
    }
    demean_mat <- function(M) {
        apply(M, 2, demean_vec)
    }

    y_dm <- demean_vec(y_raw)
    D_dm <- demean_vec(D_raw)
    X_out_dm <- demean_mat(X_out)
    X_treat_dm <- if (!is.null(X_treat)) demean_mat(X_treat) else NULL

    # Map design matrix columns to parameter roles
    # attr(model.matrix, "assign") maps each column to its term index.
    # assign == 0 is the intercept. After removing the intercept column,
    # we must also drop the 0 entry from assign.
    outcome_assign_raw <- attr(X_out_raw, "assign")
    intercept_mask <- outcome_assign_raw != 0
    outcome_assign <- outcome_assign_raw[intercept_mask]

    outcome_term_labels <- attr(
        terms(parsed$outcome_fml, data = data),
        "term.labels"
    )
    col_to_term <- outcome_term_labels[outcome_assign]

    # rho_cols: named list, one integer vector per Y-lag (ordered as lag_y)
    rho_cols <- lapply(parsed$rho_terms, function(rt) which(col_to_term == rt))
    names(rho_cols) <- lag_y

    tau_cols <- which(col_to_term %in% parsed$tau_terms)
    beta1_cols <- which(col_to_term %in% parsed$beta1_terms)

    rho2_cols <- list() # Y-lag cols in treatment equation (named by lag_y)
    lag_d_cols <- list() # D-lag cols in treatment equation (named by lag_d, placeholder)
    beta2_cols <- integer(0)

    if (parsed$has_treatment_eq && !is.null(X_treat)) {
        treat_assign_raw <- attr(X_treat_raw, "assign")
        treat_mask <- treat_assign_raw != 0
        treat_assign <- treat_assign_raw[treat_mask]

        treat_term_labels <- attr(
            terms(parsed$treatment_fml, data = data),
            "term.labels"
        )
        treat_col_to_term <- treat_term_labels[treat_assign]

        # rho2_cols: one entry per Y-lag present in treatment equation
        if (length(parsed$rho2_terms) > 0) {
            rho2_cols <- lapply(
                parsed$rho2_terms,
                function(rt) which(treat_col_to_term == rt)
            )
            names(rho2_cols) <- parsed$rho2_terms
        }

        # lag_d_cols: one entry per D-lag present in treatment equation (placeholder)
        if (length(parsed$lag_d_terms) > 0) {
            lag_d_cols <- lapply(
                parsed$lag_d_terms,
                function(lt) which(treat_col_to_term == lt)
            )
            names(lag_d_cols) <- parsed$lag_d_terms
        }

        if (length(parsed$beta2_terms) > 0) {
            beta2_cols <- which(treat_col_to_term %in% parsed$beta2_terms)
        }
    }

    list(
        y_dm = y_dm,
        D_dm = D_dm,
        X_out_dm = X_out_dm,
        X_treat_dm = X_treat_dm,
        panel_vec = panel_vec,
        time_vec = data[[time_id]],
        mu_W = mu_W,
        rho_cols = rho_cols, # list, length = n_lags_Y
        tau_cols = tau_cols,
        beta1_cols = beta1_cols,
        rho2_cols = rho2_cols, # list, length = n_rho2
        lag_d_cols = lag_d_cols, # list, length = n_lag_d (placeholder)
        beta2_cols = beta2_cols
    )
}


#' Get biased OLS estimates as starting values for GMM
#' @keywords internal
get_biased_ols <- function(dm_mats, parsed) {
    ols_outcome <- lm.fit(dm_mats$X_out_dm, dm_mats$y_dm)
    oc <- ols_outcome$coefficients

    # rho_cols is a list (one entry per Y-lag); unlist preserves order
    all_rho_cols <- unlist(dm_mats$rho_cols)

    # Assemble in θ order: [ρ₁,...,ρ_{nL}, τ terms, β₁ terms, ρ₂ terms, β₂ terms]
    theta0 <- c(
        oc[all_rho_cols],
        oc[dm_mats$tau_cols],
        oc[dm_mats$beta1_cols]
    )

    if (parsed$has_treatment_eq && !is.null(dm_mats$X_treat_dm)) {
        ols_treat <- lm.fit(dm_mats$X_treat_dm, dm_mats$D_dm)
        tc <- ols_treat$coefficients
        all_rho2_cols <- unlist(dm_mats$rho2_cols)
        # lag_d_cols placeholder: unlist will be empty when n_lag_d == 0
        all_lag_d_cols <- unlist(dm_mats$lag_d_cols)
        theta0 <- c(
            theta0,
            tc[all_rho2_cols],
            tc[all_lag_d_cols],
            tc[dm_mats$beta2_cols]
        )
    }

    names(theta0) <- build_theta_names(parsed, dm_mats)
    theta0
}


#' Build human-readable names for the θ vector
#' @keywords internal
build_theta_names <- function(parsed, dm_mats) {
    out_colnames <- colnames(dm_mats$X_out_dm)
    nms <- c(
        out_colnames[unlist(dm_mats$rho_cols)],
        out_colnames[dm_mats$tau_cols],
        out_colnames[dm_mats$beta1_cols]
    )
    if (parsed$has_treatment_eq && !is.null(dm_mats$X_treat_dm)) {
        treat_colnames <- colnames(dm_mats$X_treat_dm)
        treat_nms <- c(
            treat_colnames[unlist(dm_mats$rho2_cols)],
            treat_colnames[unlist(dm_mats$lag_d_cols)],
            treat_colnames[dm_mats$beta2_cols]
        )
        nms <- c(nms, paste0("treat_eq:", treat_nms))
    }
    nms
}
