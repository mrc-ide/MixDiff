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

n_per_group <- rep(10, n_groups)
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

# Run MCMC
# MCMC_result <- RunMCMC(obs_dat,
#                        MCMC_settings,
#                        hyperparameters,
#                        index_dates)

# Functions used within RunMCMC:
# check_MCMC_settings() - done
# initialise_aug_data() - done, but may revisit
#   -> compute_index_dates_order() - done
#   -> are_dates_incompatible() - done
# find_range() - done
# lposterior_total()
#   -> LL_total()
#         -> LL_observation_term()
#               -> LL_observation_term_by_group_delay_and_indiv()
#         -> LL_error_term()
#         -> LL_delays_term()
#               -> compute_delta()
#               -> LL_delays_term_by_group_delay_and_indiv()
#                     -> DiscrGamma()
#   -> lprior_total()
#         -> lprior_prob_error()
#         -> lprior_params_delay()
# move_Di()
#   -> propose_new_delay()
#         -> choose_delay()
#         -> discr_gamma_sample()
#   -> LL_observation_term_by_group_delay_and_indiv() - same as above
#   -> LL_error_term_by_group_delay_and_indiv()
#   -> LL_delays_term_by_group_delay_and_indiv() - same as above
#   -> get_correct_factor_new_delay()
#         -> DiscrGamma()
# move_Ei()
#   -> find_range() - same as above
#   -> propose_move_from_E0_to_E1()
#         -> find_params_gamma()
#         -> discr_gamma_sample() - same as above
#   -> compute_p_accept_move_from_E0_to_E1()
#         -> LL_observation_term_by_group_delay_and_indiv() - same as above
#         -> LL_error_term_by_group_delay_and_indiv() - same as above
#         -> LL_delays_term_by_group_delay_and_indiv() - same as above
#         -> lposterior_total() - same as above
#         -> DiscrGamma() - same as above
#         -> find_correction_factor()
#   -> propose_move_from_E1_to_E0()
#   -> compute_p_accept_move_from_E1_to_E0()
#         -> LL_observation_term_by_group_delay_and_indiv() - same as above
#         -> LL_error_term_by_group_delay_and_indiv() - same as above
#         -> LL_delays_term_by_group_delay_and_indiv() - same as above
#         -> lposterior_total() - same as above
#         -> DiscrGamma() - same as above
#         -> find_correction_factor_2() - merge with find_correction_factor()?
# find_Eis_to_swap()
# swap_Ei()
#   -> propose_move_from_E1_to_E0()
#   -> propose_move_from_E0_to_E1()
#   -> compute_p_accept_move_from_E0_to_E1()
#   -> propose_new_delay()
#   -> get_correct_factor_new_delay()
#   -> compute_delta()
#   -> LL_observation_term_by_group_delay_and_indiv()
#   -> LL_error_term_by_group_delay_and_indiv()
#   -> LL_delays_term_by_group_delay_and_indiv()
#   -> lposterior_total()
# move_zeta_gibbs()
#   -> compute_n_errors()
# move_lognormal()
#   -> lprior_params_delay()
#   -> compute_delta_group_delay_and_indiv()
#   -> LL_delays_term_by_group_delay_and_indiv()

# ----------------------------------------------------------------------------

test_that("check_MCMC_settings detects errors", {
  
  # Make burnin > n_iter
  burnin_error <- MCMC_settings
  burnin_error$chain_properties$burnin <- 501
  expect_error(check_MCMC_settings(burnin_error, index_dates),
               "Burnin must be <= n_iter")
  
  # Add an extra group to sdlog_mu
  sdlog_mu_length_error <- MCMC_settings
  sdlog_mu_length_error$moves_options$sdlog_mu[[5]] <- c(0.25, 0.25)
  expect_error(check_MCMC_settings(sdlog_mu_length_error, index_dates),
               "sdlog_mu does not have the correct length")
  
  # Add an extra group to sdlog_CV
  sdlog_CV_length_error <- MCMC_settings
  sdlog_CV_length_error$moves_options$sdlog_CV[[5]] <- c(0.25, 0.25)
  expect_error(check_MCMC_settings(sdlog_CV_length_error, index_dates),
               "sdlog_CV does not have the correct length")
  
  # Remove one sdlog_mu from one of the groups
  sdlog_mu_structure_error <- MCMC_settings
  sdlog_mu_structure_error$moves_options$sdlog_mu[[4]] <- c(0.25, 0.25)
  expect_error(check_MCMC_settings(sdlog_mu_structure_error, index_dates),
               "sdlog_mu does not have the correct structure")
  
  # Remove one sdlog_CV from one of the groups
  sdlog_CV_structure_error <- MCMC_settings
  sdlog_CV_structure_error$moves_options$sdlog_CV[[4]] <- c(0.25, 0.25)
  expect_error(check_MCMC_settings(sdlog_CV_structure_error, index_dates),
               "sdlog_CV does not have the correct structure")
  
})


# ----------------------------------------------------------------------------

# Tests for initialise_aug_dat (move to test-InitMCMC.R)

test_that("Missing dates in obs_dat are imputed during initialisation and E is given -1", {
 
  obs_dat <- list(matrix(c(
    10, 16, NA,
    12, 20, 27
  ), nrow = 2, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 15))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # No NAs in output for D
  testthat::expect_false(any(is.na(out$D[[1]])))
  
  # The missing date is given -1 in E
  testthat::expect_equal(out$E[[1]][1, 3], -1)
   
})

testthat::test_that("Dates involved in delays that are too long/short are modified", {
  # Make delays outside the bounds of mindelay = 3 and maxdelay = 10
  obs_dat <- list(matrix(c(
    10, 11, 20, # delay between 1->2 = 1 (< mindelay)
    12, 20, 31 # delay between 2->3 = 11 (> maxdelay)
  ), nrow = 2, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # Expect either [1,1] or [1,2] and [2,2] or [2,3] to be flagged as an error
  testthat::expect_true(any(out$E[[1]][1, 1] == 1 | out$E[[1]][1, 2] == 1))
  testthat::expect_true(any(out$E[[1]][2, 2] == 1 | out$E[[1]][2, 3] == 1))
})

testthat::test_that("Date involved in >1 problematic delay is modified", {
  obs_dat <- list(matrix(c(
    10, 80, 14, 18,  # date [1,2] involved in 3 problematic delays
    12, 20, 25, 27 
  ), nrow = 2, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(2, 4)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 15))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # Expect date [1,2] to be modified and flagged as error
  testthat::expect_true(out$D[[1]][1, 2] != obs_dat[[1]][1, 2])
  testthat::expect_true(out$E[[1]][1, 2] == 1)
})

testthat::test_that("Date furthest from median is corrected when only one problematic delay exists", {
  # Create only one problematic delay - delay 1->4 is 90 days (> maxdelay)
  obs_dat <- list(matrix(c(10, 15, 25, 100), nrow = 1, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(1, 4)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 15))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # Check that the date furthest from the median (100) is corrected
  # and that the other date in the problematic delay (10) is not.
  testthat::expect_true(out$D[[1]][1, 4] != obs_dat[[1]][1, 4])
  testthat::expect_true(out$D[[1]][1, 1] == obs_dat[[1]][1, 1])
  
  # Check that [1, 4] is flagged as an error and the other date [1, 1] is not
  testthat::expect_equal(out$E[[1]][1, 4], 1)
  testthat::expect_equal(out$E[[1]][1, 1], 0)
})

testthat::test_that("Transitive error is corrected", {
  # transitive delay 1->4 is 70 days (> maxdelay * 3)
  obs_dat <- list(matrix(c(10, 15, 20, 80), nrow = 1, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(3, 4)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 20))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # Check that the date furthest from the median (100) is corrected
  # and that the other date in the problematic delay (10) is not.
  testthat::expect_true(out$D[[1]][1, 4] != obs_dat[[1]][1, 4])
  
  # Check for correct error indicator
  testthat::expect_equal(out$E[[1]][1, 4], 1)
  testthat::expect_equal(out$E[[1]][1, 1], 0)
})

test_that("No dates are modified if there are no incompatible delays", {
  obs_dat <- list(matrix(c(10, 15, 20), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # Output should be the same as the input
  testthat::expect_equal(out$D[[1]], obs_dat[[1]])
  
  # E should be all 0
  testthat::expect_equal(out$E[[1]][1, ], c(0, 0, 0))
})

testthat::test_that("Two groups are handled correctly", {
  # Group 1 has an incompatible delay (1 day, mindelay = 3)
  # Group 2 has a valid delay (5 days)
  obs_dat <- list(
    matrix(c(10, 11), nrow = 1),
    matrix(c(10, 15), nrow = 1)
  )
  index_dates <- list(
    cbind(c(1, 2)),
    cbind(c(1, 2))
  )
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # Check that exactly one error is flagged in group 1
  testthat::expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 1)
  
  # Check that group 2 has no errors
  testthat::expect_true(all(out$E[[2]] == 0))
})

testthat::test_that("Mixed problematic groups are handled correctly", {
  # Group 1 has a problematic delay (> maxdelay)
  # Group 2 has a missing date
  obs_dat <- list(
    matrix(c(10, 50), nrow = 1),
    matrix(c(10, NA), nrow = 1)
  )
  index_dates <- list(
    cbind(c(1, 2)),
    cbind(c(1, 2))
  )
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # Check that exactly one error is flagged in group 1
  testthat::expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 1)
  
  # Check that missing value in group 2 is imputed and flagged as -1
  testthat::expect_false(is.na(out$D[[2]][1, 2]))
  testthat::expect_equal(out$E[[2]][1, 2], -1)
})

testthat::test_that("All delays less than mindelay handled correcty", {
  # Both specified delays are incompatible: 1->2 (delay = 1), 2->3 (delay = 2)
  # Transitive delay also incompatible: 1->3 (delay = 2)
  obs_dat <- list(matrix(c(10, 11, 12), nrow = 1, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  testthat::expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 2)
})

# Edge case to discuss:
testthat::test_that("Specified delays invalid, transitive delay valid", {
  # Both specified delays are incompatible: 1->2 (delay = 1), 2->3 (delay = 2)
  # Under previous rules the transitive delay was compatible: 1->3 (delay = 3)
  # This has been changed (now 1->3 is mindelay * 2)
  obs_dat <- list(matrix(c(10, 11, 13), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))

  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  testthat::expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 2)

  # Under previous rules date 2 is unchanged because although date 2 is
  # flagged as being in both bad pairs and set to NA, in the "missing" part
  # date 2 is constrained as both "before" date 3 and "after" date 1, so
  # floor(median(c(10, 13))) = 11, the same date as the original and so E = 0
})

testthat::test_that("All delays exceeding maxdelay handled correctly", {
  # Both specified delays exceed maxdelay: 1->2 (delay = 11), 2->3 (delay = 11)
  # Transitive delay also incompatible: 1->3 (delay = 22)
  obs_dat <- list(matrix(c(10, 21, 32), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  testthat::expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 2)
})

testthat::test_that("Missing dates in multiple groups handled correctly", {
  index_dates <- list(
    matrix(c(1, 2), nrow = 2),
    cbind(c(1, 2), c(1, 3)),
    cbind(c(1, 2), c(2, 3),c(1, 4)),
    cbind(c(1, 2), c(2, 3), c(1, 4))
  )

  obs_dat <- list(
    matrix(c(20, 24,
             19, 23), nrow = 2, byrow = TRUE),
    matrix(c(30, 35, NA,
             NA, 23, 40), nrow = 2, byrow = TRUE),
    matrix(c(17, 23, 27, NA,
             27, 33, 41, 40), nrow = 2, byrow = TRUE),
    matrix(c(23, 31, NA, 27,
             8, 10, 28, 13), nrow = 2, byrow = TRUE)
  )
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 30))

  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  expected_E <- list(
    matrix(0, nrow = 2, ncol = 2),
    matrix(c(0, 0, -1, -1, 0, 0), nrow = 2, byrow = TRUE),
    matrix(c(0, 0, 0, -1, 0, 0, 0, 0), nrow = 2, byrow = TRUE),
    matrix(c(0, 0, -1, 0, 0, 0, 0, 0), nrow = 2, byrow = TRUE)
  )

  testthat::expect_equal(out$E, expected_E)
  
})


# ----------------------------------------------------------------------------

# Compute_index_dates_order() - move to utiities tests

testthat::test_that("Simple transitivity is computed correctly", {
  # With input 1->2 and 2->3 the output should include the transitive rule 1->3
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  out <- compute_index_dates_order(index_dates)
  
  expected <- list(cbind(c(1, 2), c(2, 3), c(1, 3)))
  
  testthat::expect_equal(out[[1]], expected[[1]])
})

testthat::test_that("Multiple transitive relationships computed correctly", {
  # With input 1->2, 2->3, 3->4 output should include 1->3, 2->4, and 1->4
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(3, 4)))
  out <- compute_index_dates_order(index_dates)
  
  expected <- list(cbind(c(1, 2), c(2, 3), c(3, 4), c(1, 3), c(2, 4), c(1, 4)))
  
  testthat::expect_equal(out[[1]], expected[[1]])
})

testthat::test_that("No transitive relationships handled correctly", {
  # No transitive relationships so output should be the same as input
  index_dates <- list(cbind(c(1, 2), c(3, 4)))
  out <- compute_index_dates_order(index_dates)
  
  testthat::expect_equal(out[[1]], index_dates[[1]])
})

testthat::test_that("Duplicates in the input are removed", {
  # The duplicate 1->2 should be removed (and transitive 1->3 should be added)
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(1, 2)))
  out <- compute_index_dates_order(index_dates)
  
  expected <- list(cbind(c(1, 2), c(2, 3), c(1, 3)))
  
  testthat::expect_equal(out[[1]], expected[[1]])
})

testthat::test_that("Multiple groups are handled correctly", {
  index_dates <- list(
    matrix(c(1, 2), nrow = 2), # no change needed
    cbind(c(1, 2), c(1, 3)), # no change needed
    cbind(c(1, 2), c(2, 3), c(1, 4)), # need to add 1->3
    cbind(c(1, 2), c(2, 3), c(1, 4)) # need to add 1->3
  )
  out <- compute_index_dates_order(index_dates)
  
  expected <- list(
    matrix(c(1, 2), nrow = 2), # no change
    cbind(c(1, 2), c(1, 3)), # no change
    cbind(c(1, 2), c(2, 3), c(1, 4), c(1, 3)), # add 1->3
    cbind(c(1, 2), c(2, 3), c(1, 4), c(1, 3)) # add 1->3
  )
  
  testthat::expect_equal(out, expected)
})

# ----------------------------------------------------------------------------

# find where are_dates_incompatible function lives and move test

test_that("are_dates_incompatible works correctly", {
  mindelay <- 3
  maxdelay <- 10
  expect_true(are_dates_incompatible(10, 11, mindelay, maxdelay)) # too short
  expect_true(are_dates_incompatible(11, 13, mindelay, maxdelay)) # too short
  expect_false(are_dates_incompatible(10, 13, mindelay, maxdelay)) # mindelay
  expect_false(are_dates_incompatible(10, 20, mindelay, maxdelay)) # maxdelay
  expect_true(are_dates_incompatible(10, 21, mindelay, maxdelay)) # too long
})

# ----------------------------------------------------------------------------

# move to utilities tests

test_that("find_range finds the correct min and max across 2 groups", {
  obs_dat <- list(
    matrix(c(10, 20, 15, 25), nrow = 2),
    matrix(c(5, 50, 30, 40), nrow = 2)
  )
  
  expected_range <- c(5, 50)
  
  expect_equal(find_range(obs_dat), expected_range)
})

test_that("find_range correctly handles NA values", {
  # NAs should be ignored when other dates available
  obs_dat_na <- list(
    matrix(c(10, NA, 15, 25), nrow = 2),
    matrix(c(5, 50, NA, 40), nrow = 2)
  )
  expect_equal(find_range(obs_dat_na), c(5, 50))
  
  # a date column in one group has all NAs
  obs_dat_na_col <- list(
    matrix(c(10, 20, NA, NA), nrow = 2),
    matrix(c(5, 50, NA, 40), nrow = 2)
  )
  expect_error(
    find_range(obs_dat_na_col),
    "All dates in column 2, group 1 are NA. Check that each date column in obs_dat has >=1 date."
    )
  
  # all values are NA
  obs_dat_all_na <- list(
    matrix(c(NA, NA, NA, NA), nrow = 2)
  )
  expect_error(
    find_range(obs_dat_all_na),
    "All dates in column 1, group 1 are NA. Check that each date column in obs_dat has >=1 date."
  )
})

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
