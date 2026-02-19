g <- function(
    theta,
    data
) {
    #=====
    # resid: the moments are functions of the residuals
    #=====

    data[,
        resids_mod_1 := get(y_var) -
            theta[1] * get(lag_y_variable) -
            theta[2] * get(treatment_variable)
    ]
    data[,
        resids_mod_2 := get(treatment_variable) - theta[3] * get(lag_y_variable)
    ]

    # Calculating the standard error Sigma
    sigma_k_t <- N_T /
        (N * (N_T - 1)) *
        data[, .(resids = sum(resids_mod_1^2)), .(time)]$resids
    sigma_k_t_u <- N_T /
        (N * (N_T - 1)) *
        data[, .(resids = sum(resids_mod_2^2)), .(time)]$resids

    #=====
    # bias correction
    #=====
    # aggregate parameter saved for convenience
    phi <- (theta[1] + theta[2] * theta[3])
    # inside sum of bias correction
    l_function <- function(l) {
        sum((phi)^(0:l)) * sigma_k_t[N_T - 1 - l]
    }

    l_function_u <- function(l) {
        sum((phi)^(0:l)) * sigma_k_t_u[N_T - 1 - l]
    }
    # outside sum of bias correction
    inside_peren <- sum(mapply(l_function, 0:(N_T - 2)))
    inside_peren_u <- sum(mapply(l_function_u, 0:(N_T - 2)))

    b_1 <- -1 / N_T * inside_peren
    b_2 <- -1 / N_T * theta[3] * inside_peren
    b_3 <- -1 / N_T * theta[2] * inside_peren_u
    #=====
    # bias corrected moments
    #=====

    m1 <- data[, get(lag_y_variable)] *
        data$resids_mod_1 -
        b_1
    m2 <- data[, get(treatment_variable)] *
        data$resids_mod_1 -
        b_2
    m3 <- data[, get(lag_y_variable)] * x$resids_mod_2 #- b_3

    return(cbind(m1, m2, m3))
}
