#' @importFrom fixest demean, xpd
dbc <- function(
    fml,
    lag_y_variable,
    treatment_variable,
    time_variable,
    group_variable,
    data
) {
    fixest_fit <- fixest::feols(fml, data)

    y_var <- all.vars(fml)[1]

    demeaned_data <- fixest::demean(fixest_fit)
    demeaned_data <- as.data.table(demeaned_data)
    demeaned_data[,
        time := data[, get(time_variable)]
    ]

    N_T <- max(data[, ..time_variable])
    N <- data[, uniqueN(get(group_variable))]

    coef_est <- fixest_fit$coefficients
    x_variable_names <- names(coef_est)
    non_corrected_variables <- x_variable_names[
        !x_variable_names %in% c(lag_y_variable, treatment_variable)
    ]
    non_corrected_coefficients <- coef_est[
        !x_variable_names %in% c(lag_y_variable, treatment_variable)
    ]
    if (length(non_corrected_variables) > 0) {
        offset <- sum(non_corrected_coefficients * non_corrected_variables)
    } else {
        offset <- 0
    }

    gmm_res <- gmm::gmm(
        g,
        x = demeaned_data,
        t0 = c(0, 0, 0),
        method = "Nelder-Mead",
        control = list(reltol = 1e-25, maxit = 20000)
    )
}
library(data.table)
data <- DGP()
fml <- y ~ lag_y + D | a

lag_y_variable <- "lag_y"
treatment_variable <- "D"
time_variable <- "time"
group_variable <- "fixed"

dbc(
    fml,
    lag_y_variable,
    treatment_variable,
    time_variable,
    group_variable,
    data
)
