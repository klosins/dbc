set.seed(12345)

# Generator data with no lag Y in treatment formula
data <- DGP(N = 500, N_T = 4, rho1 = 0.2, rho2 = 0, tau = 0.5)

result <- dbc(
    outcome_fml = y ~ lag_y + D,
    # treatment_fml = D ~ lag_y,
    lag_y = "lag_y",
    treatment = "D",
    panel_id = "panel_id",
    time_id = "time",
    data = data,
    reltol = 1e-12
)

test_that("exogenous model convergence and estimates", {
    # Test the convergence of the debiasing
    # Test gmm convergence code is 0
    expect_equal(result$gmm_fit$algoInfo$convergence, 0)
    # Test that debiased estimates are close enough to true
    expect_equal(unname(abs(result$coefficients[1] - 0.2) < 1e-1), TRUE)
    expect_equal(unname(abs(result$coefficients[2] - 0.5) < 1e-1), TRUE)
    # Did we just get two parameter estimates back?
    expect_equal(length(result$coefficients), 2)
})
