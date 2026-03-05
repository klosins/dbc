set.seed(12345)

# Generator data with lag Y in treatment formula
data <- DGP(
    N = 500,
    N_T = 4,
    rho1 = 0.2,
    rho2 = 0.3,
    tau = c(0.8, 1.2),
    n_X1 = 1,
    n_X2 = 1,
    n_W = 2
)

result <- dbc(
    outcome_fml = y ~ lag_y + D:W_1 + D:W_2 + X1_1,
    treatment_fml = D ~ lag_y + X2_1,
    lag_y = "lag_y",
    treatment = "D",
    panel_id = "panel_id",
    time_id = "time",
    data = data,
    reltol = 1e-12
)

test_that("interaction model convergence and estimates", {
    # Test the convergence of the debiasing
    # Test gmm convergence code is 0
    expect_equal(result$gmm_fit$algoInfo$convergence, 0)
    # Test that debiased estimates are close enough to true
    expect_equal(unname(abs(result$coefficients[1] - 0.2) < 1e-1), TRUE)
    expect_equal(unname(abs(result$coefficients[2] - 0.8) < 1e-1), TRUE)
    expect_equal(unname(abs(result$coefficients[3] - 1.2) < 1e-1), TRUE)
    expect_equal(unname(abs(result$coefficients[4] - 1) < 1e-1), TRUE)
    expect_equal(unname(abs(result$coefficients[5] - 0.3) < 1e-1), TRUE)
    expect_equal(unname(abs(result$coefficients[6] - 1) < 1e-1), TRUE)
    # Did we just get three parameter estimates back?
    expect_equal(length(result$coefficients), 6)
})
