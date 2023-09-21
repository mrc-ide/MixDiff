test_that("RunMCMC checks for convergence", {

  obs <- simulate_single_group()
  settings <- bad_settings_one_group()
  out <- RunMCMC(
    obs$obs_dat,
    settings$mcmc_settings,
    settings$hyperparameters,
    obs$index_dates
  )
  expect_named(out$convergence, c("zeta", "mu", "CV"))

  expect_vector(out$convergence$zeta, size = 1)
  expect_vector(out$convergence$mu, size = 1)
  expect_vector(out$convergence$CV, size = 1)

  expect_vector(out$convergence$mu[[1]], size = 1)
  expect_vector(out$convergence$CV[[1]], size = 1)


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
  ## Here we can expect convergence
  expect_mapequal(out$convergence, ptype)
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
  ## Check for structure here because one of the flags might be TRUE
  ## just by chance.
  expect_named(out$convergence, c("zeta", "mu", "CV"))

  expect_vector(out$convergence$zeta, size = 1)
  expect_vector(out$convergence$mu, size = 4)
  expect_vector(out$convergence$CV, size = 4)

  expect_vector(out$convergence$mu[[1]], size = 1)
  expect_vector(out$convergence$mu[[2]], size = 2)
  expect_vector(out$convergence$mu[[3]], size = 3)
  expect_vector(out$convergence$mu[[3]], size = 3)

  expect_vector(out$convergence$CV[[1]], size = 1)
  expect_vector(out$convergence$CV[[2]], size = 2)
  expect_vector(out$convergence$CV[[3]], size = 3)
  expect_vector(out$convergence$CV[[3]], size = 3)

})
