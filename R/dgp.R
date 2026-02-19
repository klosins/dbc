#' Generate data from simple dynamic panel
#'
#' @param N Number of individuals
#' @param N_T Number of time periods
#' @param rho_1 Coefficient on lagged outcome in outcome equation
#' @param rho_2 Coefficient on lagged outcome in treatment equation
DGP <- function(N = 100, N_T = 10, rho_1 = .8, rho_2 = 0, tau = 1) {
    # Drawing fixed effects
    a <- rnorm(N, 0, 5)

    # ID number for person
    fixed <- rep(1:N, each = N_T) # fixed effect dummies

    # Time label
    time <- rep(1:N_T, N)

    # Creating base data
    DT <- data.table(fixed, time, a = rep(a, each = N_T))

    # base treatment that we will write over to make endogenous
    DT$D <- rnorm(N * N_T, 0, 1)
    #generate start y
    DT[time == 1, y := a + tau * D]

    # creating lag structure
    for (year in c(2:N_T)) {
        DT[, lag_y := shift(y), .(fixed)]
        DT[time == year, D := rho_2 * lag_y + rnorm(.N, 0, 1)]
        DT[time == year, y := a + tau * D + rho_1 * lag_y + rnorm(.N, 0, 1)]
    }
    #format data
    DT <- DT[time > 1, ] # one lag removes one time period
    # Performing the within transformation on all columns except 'id' and 'time'
    # DT_within <- copy(DT)
    # cols_to_transform <- setdiff(names(DT), c("fixed", "time"))
    # DT_within[,
    #     (cols_to_transform) := lapply(.SD, function(col) col - mean(col)),
    #     by = fixed,
    #     .SDcols = cols_to_transform
    # ]
    # returning the data objects
    return(DT)
}
