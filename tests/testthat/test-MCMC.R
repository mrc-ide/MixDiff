test_that("RunMCMC checks for convergence", {

  obs <- simulate_single_group()
  settings <- bad_settings_one_group()
  out <- RunMCMC(
    obs,
    settings$mcmc_settings,
    settings$hyperparameters,
    index_dates
  )
  ptype <- list(zeta = NULL, mu = NULL, CV = NULL)
  expect_vector(
    out$convergence, ptype = ptype
  )
})
