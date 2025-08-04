# # Load simulated data ---------------------------------------------------------
# sim_data <- readRDS("tests/testdata/sim_data.rds")
# params <- readRDS("tests/testdata/params.rds")

# Simulate data to use
n_groups <- 4
index_dates <- list(
  matrix(c(1, 2), nrow = 2),
  cbind(c(1, 2), c(1, 3)),
  cbind(c(1, 2), c(2, 3),c(1, 4)),
  cbind(c(1, 2), c(2, 3), c(1, 4))
)
theta <- list(
  mu = list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13)),
  CV = list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5)),
  prop_missing_data = 0.2,
  zeta = 0.05
)

n_per_group <- rep(100, n_groups)
range_dates <- c(0, 30)

simul_dat <- simul_true_data(theta, n_per_group, range_dates, index_dates,
                             simul_error = TRUE)
obs_dat <- simul_dat$obs_dat

# Set up hyperparameters
hyperparameters <- list(
  shape1_prob_error = 3,
  shape2_prob_error = 12,
  mean_mean_delay = 10,
  mean_CV_delay = 10)

# Set up MCMC
MCMC_settings <- list(
  moves_switch = list(D_on = TRUE, E_on = TRUE, swapE_on = TRUE,
                      mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
  moves_options = list(
    fraction_Di_to_update = 1 / 10,
    move_D_by_groups_of_size = 1,
    fraction_Ei_to_update = 1 / 10,
    sdlog_mu = list(
      0.05,
      c(0.15, 0.15),
      c(0.15, 0.15, 0.15),
      c(0.25, 0.25, 0.25)
    ),
    sdlog_CV = list(
      0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
  ),
  init_options = list(
    mindelay = 0,
    maxdelay = 20
  ),
  chain_properties = list(
    n_iter = 500,
    burnin = 50,
    record_every = 10
  )
)

# Initialise augmented data ---------------------------------------------------

# obs_dat <- sim_data$obs_dat
# index_dates <- params$index_dates
# MCMC_settings <- list(init_options = list(mindelay = 0, maxdelay = 100))

test_that("initialise_aug_data returns correctly structured output", {
  
  aug_dat <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  expect_type(aug_dat, "list")
  expect_named(aug_dat, c("D", "E"))
  expect_length(aug_dat$D, length(sim_data$obs_dat))
  expect_equal(dim(aug_dat$D[[1]]), dim(sim_data$obs_dat[[1]]))
  
})


