test_that("RunMCMC checks for convergence", {

  obs <- simulate_single_group()
  settings <- bad_settings_one_group()
  out <- RunMCMC(
    obs$obs_dat,
    settings$mcmc_settings,
    settings$hyperparameters,
    obs$index_dates
  )
  ptype <- list(zeta = NULL, mu = NULL, CV = NULL)
  expect_vector(
    out$convergence, ptype = ptype
  )
})
