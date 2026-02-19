# Feb 13th 2024
# using the R gmm package
# createing the moments for
# the treatment function of the
# past story

library(gmm)
library(data.table)
library(parallel)
library(xtable)
library(MASS)

#================
# Section 1: Set up DGP
#================

# create one dataset,
DGP <- function(N, N_T = 10, rho = .8, tau = 1) {
  rho_2 <- .6
  # drawing fixed effects
  a <- rnorm(N, 0, 5)
  # id number for person
  fixed <- rep(1:N, each = N_T) # fixed effect dummies
  # time label
  time <- rep(1:N_T, N) # time
  # creating base data
  DT <- data.table(fixed, time, a = rep(a, each = N_T))
  # base treatment that we will write over to make endogenous
  DT$D <- rnorm(N * N_T, 0, 1)
  #generate start y
  DT[time == 1, y := a + tau * D]

  # creating lag structure
  for (year in c(2:N_T)) {
    DT[, lag_y := shift(y), .(fixed)]
    DT[time == year, D := rho_2 * lag_y + rnorm(.N, 0, 1)]
    DT[time == year, y := a + tau * D + rho * lag_y + rnorm(.N, 0, 1)]
  }
  #format data
  DT <- DT[time > 1, ] # one lag removes one time period
  # Performing the within transformation on all columns except 'id' and 'time'
  DT_within <- copy(DT)
  cols_to_transform <- setdiff(names(DT), c("fixed", "time"))
  DT_within[,
    (cols_to_transform) := lapply(.SD, function(col) col - mean(col)),
    by = fixed,
    .SDcols = cols_to_transform
  ]
  # returning the data objects
  return(DT_within)
}

#================
# Section 2: Set Up Moments
#================

# theta = (rho_1, tau , rho_2)

# Define moment conditions
g <- function(theta, DT) {
  # save parameters from data
  N_T <- max(DT$time)
  N <- max(DT$fixed)

  #=====
  # resid: the moments are functions of the residuals
  #=====
  DT[, resids_mod_1 := y - theta[1] * lag_y - theta[2] * D]
  DT[, resids_mod_2 := D - theta[3] * lag_y]

  # Calculating the standard error Sigma
  sigma_k_t <- N_T /
    (N * (N_T - 1)) *
    DT[, .(resids = sum(resids_mod_1^2)), .(time)]$resids
  sigma_k_t_u <- N_T /
    (N * (N_T - 1)) *
    DT[, .(resids = sum(resids_mod_2^2)), .(time)]$resids

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

  m1 <- DT$lag_y * DT$resids_mod_1 - b_1
  m2 <- DT$D * DT$resids_mod_1 - b_2
  m3 <- DT$lag_y * DT$resids_mod_2 #- b_3

  return(cbind(m1, m2, m3))
}


#================
# Section 3: Run GMM
#================

# Initial parameter values
theta_init <- c(0, 0, 0)

# Estimating the model
DT <- DGP(N = 100, N_T = 10, rho = .7, tau = 1)
gmm_res <- gmm(
  g,
  x = DT,
  t0 = theta_init,
  crit = 1e-25,
  method = "Nelder-Mead",
  control = list(reltol = 1e-25, maxit = 20000)
)

# Display results
summary(gmm_res)$coefficients[1:3]

#lm(D ~ lag_y, data = DT)
lm(y ~ lag_y + D, data = DT)


# Define moment conditions
g <- function(theta, DT) {
  #=====
  # resid: the moments are functions of the residuals
  #=====
  resids_mod_1 <- DT$y - theta[1] - theta[1] * DT$lag_y - theta[2] * DT$D

  #=====
  # bias correction
  #=====

  #=====
  # bias corrected moments
  #=====

  m1 <- (DT$lag_y * resids_mod_1)
  m2 <- (DT$D * resids_mod_1)

  return(cbind(m1, m2))
}
