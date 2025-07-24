# Simulate data for tests

# Parameters
n_groups <- 4
n_per_group <- rep(100, n_groups)
n_dates <- c(2, 3, 4, 4)

mu_baseline <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
cv_baseline <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))

theta <- list(
  prop_missing_data = 0.2,
  zeta = 0.05,
  mu = mu_baseline,
  CV = cv_baseline
)

range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
                             as.Date("01/01/2015", "%d/%m/%Y")))

index_dates <- list(
  matrix(c(1, 2), nrow = 2),
  cbind(c(1, 2), c(1, 3)),
  cbind(c(1, 2), c(2, 3),c(1, 4)),
  cbind(c(1, 2), c(2, 3), c(1, 4))
)

# Save parameters as a list
params <- list(
  n_groups = n_groups,
  n_per_group = n_per_group,
  n_dates = n_dates,
  theta = theta,
  range_dates = range_dates,
  index_dates = index_dates
)

saveRDS(params, "tests/testdata/params.rds")


# Simulate data
sim_data <- simul_true_data(theta,
                            n_per_group,
                            range_dates,
                            index_dates,
                            simul_error = TRUE,
                            remove_allNA_indiv = TRUE)

saveRDS(sim_data, "tests/testdata/sim_data.rds")

obs_dat <- sim_data$obs_dat

init_aug <- initialise_aug_data(obs_dat, index_dates, MCMC_settings = mcmc_settings)

theta <- initialise_theta_from_aug_dat(init_aug, index_dates)




# MCMC_settings
mcmc_settings <- list(
  # moves_switch: booleans stating whether each parameter/augmented data
  # should be moved in the procedure or not.
  moves_switch = list(
    D_on = TRUE, # augmented dates (latent true dates)
    E_on = TRUE, # error indicators (-1, 0, 1)
    swapE_on = TRUE, # swaps error indicators - explore alternative errors
    mu_on = TRUE, # mean of each delay distribution
    CV_on = TRUE, # cv of each delay
    zeta_on = TRUE # probability of error
  ),
  # moves_options: 
  moves_options = list(
    # Fraction of augmented dates to be updated at each iteration of the MCMC.
    fraction_Di_to_update = 1 / 10,
    # Number of augmented dates to be updated simultaneously in each group.
    move_D_by_groups_of_size = 1,
    # Fraction of indicators of whether observed dates are erroneous to be
    # updated at each iteration of the MCMC.
    fraction_Ei_to_update = 1 / 10,
    # List of SDs used for proposing moves of the mean delays of length n_groups.
    # Each element in the list should be a vector with length given by the
    # numbers of delays to be considered in this group.
    sdlog_mu = list(
      0.05,
      c(0.15, 0.15),
      c(0.15, 0.15, 0.15),
      c(0.25, 0.25, 0.25)
    ),
    # Same as above but for proposing moves of the CV of delays.
    sdlog_CV = list(
      0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
  ),
  # minimum and maximum delays, below/above which dates are considered
  # incompatible with one another at the initialisation stage of the MCMC.
  init_options = list(mindelay = 0, maxdelay = 100),
  # total number of iterations, initial burnin and then after burnin how many
  # iterations should be recorded (thinning). (500 - 50) / 10 = 45 samples from
  # the posterior for each dataset.
  chain_properties = list(n_iter = 500, burnin = 50, record_every = 10)
)

# Hyperparameters
hyperparameters <- list(
  # scalars giving the 1st and 2nd shape parameters for the beta prior for zeta
  shape1_prob_error = 3,
  shape2_prob_error = 12,
  # scalars giving the mean of the exponential prior used for mu and CV
  mean_mean_delay = 100,
  mean_CV_delay = 100
)




