# broom methods for "dbc" objects
#
# These methods enable compatibility with broom, modelsummary, and other
# tidy-workflow packages. broom is listed in Suggests, so these methods are
# only active when broom is loaded (lazy S3 registration via NAMESPACE).

#' Tidy a dbc object into a coefficient table
#'
#' Returns one row per parameter with columns `term`, `estimate`, `std.error`,
#' `statistic`, `p.value`, and `group`. The `group` column labels each row
#' as `"Outcome"` or `"Treatment"`, enabling vertically stacked equation
#' subpanels in `modelsummary::modelsummary()` via
#' `shape = group + term ~ model`.
#'
#' @param x A `dbc` object.
#' @param conf.int Logical: include confidence interval columns? Default FALSE.
#' @param conf.level Confidence level for intervals. Default 0.95.
#' @param ... Currently unused.
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
#' broom::tidy(fit)
#' broom::tidy(fit, conf.int = TRUE, conf.level = 0.95)
#' @exportS3Method broom::tidy
tidy.dbc <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
    coefs <- x$coefficients
    ses <- x$se
    zvals <- coefs / ses
    pvals <- 2 * stats::pnorm(abs(zvals), lower.tail = FALSE)
    terms <- names(coefs)

    # Identify treatment-equation rows and strip the prefix
    is_treat <- startsWith(terms, "treat_eq:")
    terms[is_treat] <- sub("^treat_eq:", "", terms[is_treat])

    out <- data.frame(
        group = ifelse(is_treat, "Treatment", "Outcome"),
        term = terms,
        estimate = unname(coefs),
        std.error = unname(ses),
        statistic = unname(zvals),
        p.value = unname(pvals),
        stringsAsFactors = FALSE
    )

    if (conf.int) {
        alpha <- 1 - conf.level
        z <- stats::qnorm(1 - alpha / 2)
        out$conf.low <- unname(coefs) - z * unname(ses)
        out$conf.high <- unname(coefs) + z * unname(ses)
    }

    out
}


#' Glance at a dbc object
#'
#' Returns a one-row data frame of model-level summary statistics. Compatible
#' with [broom::glance()] and `modelsummary::modelsummary()`.
#'
#' `fe_panel` is always `"Yes"` because the DBC estimator absorbs individual
#' fixed effects via within-demeaning. `fe_time` is `"No"` because time fixed
#' effects are not automatically absorbed; include time dummies as covariates
#' if needed.
#'
#' @param x A `dbc` object.
#' @param ... Currently unused.
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
#' broom::glance(fit)
#' @exportS3Method broom::glance
glance.dbc <- function(x, ...) {
    convergence <- tryCatch(
        x$gmm_fit$algoInfo$convergence == 0,
        error = function(e) NA
    )
    data.frame(
        N = x$N,
        N_T = x$N_T,
        nobs = x$N * x$N_T,
        phi = signif(x$phi, 3),
        converged = convergence,
        stringsAsFactors = FALSE
    )
}
