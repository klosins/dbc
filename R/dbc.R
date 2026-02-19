#' @importFrom fixest demean, xpd
dbc <- function(
    fml,
    lag_y_variable,
    treatment_variable,
    time_variable,
    group_variable,
    data
) {
    # First, use fixest to fit the biased model
    fixest_fit <- fixest::feols(fml, data)

    # Get the outcome variable
    y_var <- all.vars(fml)[1]

    # Do the demeaning, and add back on the time variable
    demeaned_data <- fixest::demean(fixest_fit)
    demeaned_data <- as.data.table(demeaned_data)
    demeaned_data[,
        time := data[, get(time_variable)]
    ]

    # Get group/time stats
    N_T <- max(data[, ..time_variable])
    N <- data[, uniqueN(get(group_variable))]

    # Extract out estimates from the fixest package
    # To my mind, this would be most useful for a model with lots of control X variables
    coef_est <- fixest_fit$coefficients
    x_variable_names <- names(coef_est)

    # Extract out these control variables and find their XB vector values
    non_corrected_variables <- x_variable_names[
        !x_variable_names %in% c(lag_y_variable, treatment_variable)
    ]
    non_corrected_coefficients <- coef_est[
        !x_variable_names %in% c(lag_y_variable, treatment_variable)
    ]

    if (length(non_corrected_variables) > 0) {
        # X'B (but i don't think this line works at the moment)
        offset <- non_corrected_coefficients %*% data[, non_corrected_variables]
    } else {
        offset <- 0
    }

    # Then, pass this demeaned data to the GMM function
    # TODO: Hard to pass additional arguments into g through gmm, so we'd need to setup and an env with those additional arguments and environmental variables
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
