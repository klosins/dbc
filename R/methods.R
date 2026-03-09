# Coefficient table helper
.coef_table <- function(x) {
    coefs <- x$coefficients
    ses <- x$se
    zvals <- coefs / ses
    pvals <- 2 * pnorm(abs(zvals), lower.tail = FALSE)
    cbind(
        Estimate = coefs,
        `Std. Error` = ses,
        `z value` = zvals,
        `Pr(>|z|)` = pvals
    )
}

#' Print a dbc object
#'
#' @param x A `dbc` object.
#' @param digits Number of significant digits.
#' @param ... Further arguments passed to [printCoefmat()].
#' @method print dbc
#' @examples
#' set.seed(1)
#' dat <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y = "lag_y", treatment = "D",
#'     panel_id = "panel_id", time_id = "time",
#'     data = dat
#' )
#' print(fit)
#' @export
print.dbc <- function(x, digits = 4, ...) {
    cat("Dynamic Bias-Corrected Estimator\n\n")
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")

    coef_mat <- .coef_table(x)

    # Split outcome and treatment equation rows for readability
    all_nms <- rownames(coef_mat)
    treat_idx <- startsWith(all_nms, "treat_eq:")

    if (any(!treat_idx)) {
        cat("Outcome equation:\n")
        printCoefmat(
            coef_mat[!treat_idx, , drop = FALSE],
            digits = digits,
            signif.stars = TRUE,
            ...
        )
    }

    if (any(treat_idx)) {
        cat("\nTreatment equation:\n")
        rownames(coef_mat)[treat_idx] <- sub(
            "^treat_eq:",
            "",
            all_nms[treat_idx]
        )
        printCoefmat(
            coef_mat[treat_idx, , drop = FALSE],
            digits = digits,
            signif.stars = TRUE,
            ...
        )
    }

    cat(
        "\nN =",
        x$N,
        " T =",
        x$N_T,
        " phi =",
        round(x$phi, digits),
        "\n"
    )
    invisible(x)
}

#' Summary of a dbc object
#'
#' Returns a `summary.dbc` object containing the coefficient table, model
#' dimensions, and GMM convergence status.
#'
#' @param object A `dbc` object.
#' @param ... Currently unused.
#' @method summary dbc
#' @examples
#' set.seed(1)
#' dat <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y = "lag_y", treatment = "D",
#'     panel_id = "panel_id", time_id = "time",
#'     data = dat
#' )
#' summary(fit)
#' @export
summary.dbc <- function(object, ...) {
    convergence <- tryCatch(
        object$gmm_fit$algoInfo$convergence,
        error = function(e) NA_integer_
    )
    structure(
        list(
            call = object$call,
            coefficients = .coef_table(object),
            phi = object$phi,
            N = object$N,
            N_T = object$N_T,
            convergence = convergence
        ),
        class = "summary.dbc"
    )
}

#' Print a summary.dbc object
#'
#' @param x A `summary.dbc` object.
#' @param digits Number of significant digits.
#' @param ... Further arguments passed to [printCoefmat()].
#' @method print summary.dbc
#' @examples
#' set.seed(1)
#' dat <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y = "lag_y", treatment = "D",
#'     panel_id = "panel_id", time_id = "time",
#'     data = dat
#' )
#' print(summary(fit))
#' @export
print.summary.dbc <- function(x, digits = 4, ...) {
    cat("Dynamic Bias-Corrected Estimator\n\n")
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")

    all_nms <- rownames(x$coefficients)
    treat_idx <- startsWith(all_nms, "treat_eq:")

    if (any(!treat_idx)) {
        cat("Outcome equation:\n")
        printCoefmat(
            x$coefficients[!treat_idx, , drop = FALSE],
            digits = digits,
            signif.stars = TRUE,
            ...
        )
    }

    if (any(treat_idx)) {
        cat("\nTreatment equation:\n")
        rownames(x$coefficients)[treat_idx] <- sub(
            "^treat_eq:",
            "",
            all_nms[treat_idx]
        )
        printCoefmat(
            x$coefficients[treat_idx, , drop = FALSE],
            digits = digits,
            signif.stars = TRUE,
            ...
        )
    }

    cat(
        "\nN =",
        x$N,
        " T =",
        x$N_T,
        " phi =",
        round(x$phi, digits)
    )
    if (!is.na(x$convergence)) {
        conv_msg <- if (x$convergence == 0) "(converged)" else "(NOT converged)"
        cat("  GMM:", conv_msg)
    }
    cat("\n")
    invisible(x)
}


# coef, vcov
# ----------

#' Extract coefficients from a dbc object
#'
#' @param object A `dbc` object.
#' @param ... Currently unused.
#' @method coef dbc
#' @examples
#' set.seed(1)
#' dat <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y = "lag_y", treatment = "D",
#'     panel_id = "panel_id", time_id = "time",
#'     data = dat
#' )
#' coef(fit)
#' @export
coef.dbc <- function(object, ...) {
    object$coefficients
}

#' Extract the variance-covariance matrix from a dbc object
#'
#' @param object A `dbc` object.
#' @param ... Currently unused.
#' @method vcov dbc
#' @examples
#' set.seed(1)
#' dat <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y = "lag_y", treatment = "D",
#'     panel_id = "panel_id", time_id = "time",
#'     data = dat
#' )
#' vcov(fit)
#' @export
vcov.dbc <- function(object, ...) {
    object$vcov
}


# confint.dbc
# -----------

#' Confidence intervals for a dbc object
#'
#' Computes Wald-type confidence intervals using the asymptotic normal
#' distribution.
#'
#' @param object A `dbc` object.
#' @param parm  Parameters to include: a character vector of names or an
#'   integer vector of positions. Defaults to all parameters.
#' @param level Confidence level (default 0.95).
#' @param ... Currently unused.
#' @method confint dbc
#' @examples
#' set.seed(1)
#' dat <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y = "lag_y", treatment = "D",
#'     panel_id = "panel_id", time_id = "time",
#'     data = dat
#' )
#' confint(fit)
#' confint(fit, level = 0.90)
#' confint(fit, parm = "D")
#' @export
confint.dbc <- function(object, parm, level = 0.95, ...) {
    cf <- coef(object)
    ses <- object$se

    if (!missing(parm)) {
        cf <- cf[parm]
        ses <- ses[parm]
    }

    alpha <- 1 - level
    z <- qnorm(1 - alpha / 2)

    ci <- cbind(cf - z * ses, cf + z * ses)
    colnames(ci) <- paste0(
        format(100 * c(alpha / 2, 1 - alpha / 2), trim = TRUE),
        "%"
    )
    ci
}


# nobs.dbc
# --------

#' Number of observations used in a dbc model
#'
#' Returns the total number of panel observations (N x T).
#'
#' @param object A `dbc` object.
#' @param ... Currently unused.
#' @method nobs dbc
#' @examples
#' set.seed(1)
#' dat <- DGP(N = 100, N_T = 4, rho1 = 0.3, rho2 = 0.2, tau = 0.5)
#' fit <- dbc(
#'     outcome_fml   = y ~ lag_y + D,
#'     treatment_fml = D ~ lag_y,
#'     lag_y = "lag_y", treatment = "D",
#'     panel_id = "panel_id", time_id = "time",
#'     data = dat
#' )
#' nobs(fit)
#' @export
nobs.dbc <- function(object, ...) {
    object$N * object$N_T
}
