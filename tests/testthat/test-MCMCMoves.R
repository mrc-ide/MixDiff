
# Load simulated data
sim_data <- readRDS("tests/testdata/sim_data.rds")
params <- readRDS("tests/testdata/params.rds")

# Other params
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

### Move_Di

i <- 1
group_idx <- 1
date_idx <- 1
curr_aug_dat <- initialise_aug_data(sim_data$obs_dat,
                                    compute_index_dates_order(index_dates),
                                    MCMC_settings = mcmc_settings)
theta <- params$theta
obs_dat <- sim_data$obs_dat
hyperparameters <- hyperparameters
index_dates <- params$index_dates
range_dates <- find_range(obs_dat)


test <- move_Di(i,
                group_idx,
                date_idx,
                curr_aug_dat,
                theta,
                obs_dat,
                hyperparameters,
                index_dates,
                range_dates = NULL)


if (test$accept == 0) {
  expect_equal(test$new_aug_dat$D[[1]][1,1], curr_aug_dat$D[[1]][1,1])
} else {
  expect_failure(expect_equal(test$new_aug_dat$D[[1]][1,1], curr_aug_dat$D[[1]][1,1]))
}
 # new value: test$new_aug_dat$D[[1]][1,1]
 # old value: curr_aug_dat$D[[1]][1,1]


#------------------------------------------------------------------------------

test_that("xyz", {
  
})

#------------------------------------------------------------------------------

test_that("ratio_post is the same as ratio_post_long in move_Di", {
  
})

#------------------------------------------------------------------------------

test_that("ratio_post is the same as ratio_post_long in compute_p_accept_move_from_E0_to_E1", {
  
})

#------------------------------------------------------------------------------

test_that("ratio_post is the same as ratio_post_long in compute_p_accept_move_from_E1_to_E0", {
  
})

#------------------------------------------------------------------------------

test_that("ratio_post is the same as ratio_post_long in swap_Ei", {
  
})

#------------------------------------------------------------------------------

test_that("ratio_post is the same as ratio_post_long in move_lognormal", {
  
})