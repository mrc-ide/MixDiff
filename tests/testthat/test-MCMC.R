test_that("RunMCMC checks for convergence", {

  obs <- simulate_single_group()
  settings <- bad_settings_one_group()
  out <- RunMCMC(
    obs$obs_dat,
    settings$mcmc_settings,
    settings$hyperparameters,
    obs$index_dates
  )
  ptype <- list(
      zeta = FALSE,
      mu =list(list(FALSE)),
      CV = list(list(FALSE))
  )
  expect_mapequal(
    out$convergence, expected = ptype
  )

})


test_that("RunMCMC checks for convergence", {
  skip_on_cran()
  obs <- simulate_single_group()
  settings <- good_settings()
  out <- RunMCMC(
    obs$obs_dat,
    settings$mcmc_settings,
    settings$hyperparameters,
    obs$index_dates
  )
  ptype <- list(
      zeta = TRUE,
      mu =list(list(TRUE)),
      CV = list(list(TRUE))
  )
  expect_mapequal(
    out$convergence, expected = ptype
  )

})


test_that(
    "RunMCMC checks for convergence for multiple groups", {

  obs <- simulate_four_groups()
  settings <- bad_settings_four_groups()
  out <- RunMCMC(
    obs$obs_dat,
    settings$mcmc_settings,
    settings$hyperparameters,
    obs$index_dates
  )
  ptype <- list(
      zeta = FALSE,
      mu =list(list(FALSE)),
      CV = list(list(FALSE))
  )
  expect_mapequal(
    out$convergence, expected = ptype
  )

})
