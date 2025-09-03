# Functions in Utilities.R
# - DiscrGamma (x)
# - discr_gamma_sample (x)
# - date_to_int (x)
# - int_to_date (x)
# - find_params_beta (x)
# - find_params_gamma (x)
# - compute_delta_group_delay_and_indiv (x)
# - compute_delta (x)
# - find_range (x)
# - compute_index_dates_order (x)
# - check_MCMC_settings (x)

# ----------------------------------------------------------------------------

# DiscrGamma -----------------------------------------------------------------

test_that("DiscrGamma works for a known, simple case correctly", {
  # expect a non-zero probability for k=5 when mu=5
  prob <- DiscrGamma(k = 5, mu = 5, cv = 0.5, log = FALSE)
  expect_true(prob > 0 & prob < 1)
  log_prob <- DiscrGamma(k = 5, mu = 5, cv = 0.5, log = TRUE)
  expect_equal(log_prob, log(prob))
})

test_that("DiscrGamma works with vector input for k", {
  k_vec <- c(1, 5, 10)
  mu <- 5
  cv <- 0.5
  probs <- DiscrGamma(k = k_vec, mu = mu, cv = cv, log = FALSE)
  expect_length(probs, 3)
  expect_equal(probs[1], DiscrGamma(k = 1, mu = mu, cv = cv, log = FALSE))
  expect_equal(probs[2], DiscrGamma(k = 5, mu = mu, cv = cv, log = FALSE))
  expect_equal(probs[3], DiscrGamma(k = 10, mu = mu, cv = cv, log = FALSE))
})

test_that("Parameterisation with cv or sigma is equivalent", {
  mu <- 10
  cv <- 0.3
  sigma <- mu * cv
  prob_cv <- DiscrGamma(k = 8, mu = mu, cv = cv, log = FALSE)
  prob_sigma <- DiscrGamma(k = 8, mu = mu, sigma = sigma, log = FALSE)
  expect_equal(prob_cv, prob_sigma)
})

test_that("DiscrGamma handles k far from the mean", {
  prob_far <- DiscrGamma(k = 50, mu = 5, cv = 0.5, log = FALSE)
  # probability should be very close to zero
  expect_lt(prob_far, 1e-6)
})

test_that("DiscrGamma returns error messages with invalid inputs", {
  expect_error(DiscrGamma(k = 5, mu = 10, cv = -1), "cv must be >= 0.")
  expect_error(DiscrGamma(k = 5, mu = 10, sigma = -1), "sigma must be >= 0.")
})

# discr_gamma_sample ---------------------------------------------------------

test_that("discr_gamma_sample works as expected", {
  
  set.seed(10)
  delay_sim <- discr_gamma_sample(1000, mu = 5, cv = 0.5)
  expect_equal(mean(delay_sim), 5, tolerance = 0.05)
  
  delay_sim <- discr_gamma_sample(1000, mu = 7, cv = 0.5)
  expect_equal(mean(delay_sim), 7, tolerance = 0.05)
  
  delay_sim <- discr_gamma_sample(1000, mu = 10, cv = 0.5)
  expect_equal(mean(delay_sim), 10, tolerance = 0.05)
  
  delay_sim <- discr_gamma_sample(1000, mu = 12, cv = 0.5)
  expect_equal(mean(delay_sim), 12, tolerance = 0.05)
  
  delay_sim <- discr_gamma_sample(1000, mu = 15, cv = 0.5)
  expect_equal(mean(delay_sim), 15, tolerance = 0.05)
  
})

# date_to_int ----------------------------------------------------------------

test_that("date_to_int works for the origin date", {
  expect_equal(date_to_int(as.Date("1970-01-01")), 0)
})

test_that("date_to_int handles days after the origin correctly", {
  expect_equal(date_to_int(as.Date("1970-01-11")), 10)
})

# avoiding this with range_dates but just to check...
test_that("date_to_int handles days before the origin correctly", {
  expect_equal(date_to_int(as.Date("1969-12-31")), -1)
})

test_that("date_to_int works with an origin that is not the default", {
  new_origin <- as.Date("2025-01-01")
  expect_equal(date_to_int(as.Date("2025-01-11"), origin = new_origin), 10)
  expect_equal(date_to_int(as.Date("2024-12-22"), origin = new_origin), -10)
})

test_that("date_to_int returns errors for invalid input types", {
  expect_error(date_to_int("2025-01-01"),
               "`date` argument must be a Date object. Use as.Date().")
  expect_error(date_to_int(as.Date("2025-01-01"), origin = "2024-01-01"),
               "`origin` argument must be a Date object. Use as.Date().")
})

# int_to_date ----------------------------------------------------------------

test_that("int_to_date works for zero", {
  expect_equal(int_to_date(0), as.Date("1970-01-01"))
})

test_that("int_to_date works when using default origin", {
  expect_equal(int_to_date(10), as.Date("1970-01-11"))
})

test_that("int_to_date works with an origin that is not the default", {
  new_origin <- as.Date("2020-02-25")
  expect_equal(int_to_date(5, origin = new_origin), as.Date("2020-03-01"))
})

test_that("int_to_date returns errors for invalid input types", {
  expect_error(int_to_date("xyz"), "`int` argument must be numeric.")
  expect_error(int_to_date(10, origin = "2024-01-01"),
               "`origin` argument must be a Date object. Use as.Date().")
})

# find_params_beta -----------------------------------------------------------

test_that("find_params_beta correctly recovers known shape parameters", {
  known_shape1 <- 3
  known_shape2 <- 12
  # mean = shape1 / (shape1 + shape2) = 3 / 15 = 0.2
  expected_mean <- known_shape1 / (known_shape1 + known_shape2)
  # var = shape1 * shape2 / ((shape1 + shape2)^2 * (shape1 + shape2 + 1)) =
  # 36 / (15^2 * 16) = 0.01
  expected_var <- (known_shape1 * known_shape2) / 
    ((known_shape1 + known_shape2)^2 * (known_shape1 + known_shape2 + 1))
  calculated_params <- find_params_beta(mean = expected_mean, var = expected_var)
  expect_equal(calculated_params, c(known_shape1, known_shape2))
})

# find_params_gamma ----------------------------------------------------------

test_that("find_params_gamma correctly recovers known shape and scale parameters", {
  known_shape <- 4
  known_scale <- 2
  expected_mean <- known_shape * known_scale # 4 * 2 = 8
  expected_sigma <- sqrt(known_shape * known_scale^2) # sqrt(4 * 2^2) = sqrt(16) = 4
  calculated_params <- find_params_gamma(mean = expected_mean, sigma = expected_sigma)
  expect_equal(calculated_params, c(known_shape, known_scale))
})

# compute_delta_group_delay_and_indiv ----------------------------------------

D <- list(matrix(c(10, 20,
                   12, 25,
                   15, 28),
                 nrow = 3, ncol = 2, byrow = TRUE),
          matrix(c(100, 105, 108, 110,
                   110, 120, 125, 130,
                   130, 145, 150, 140),
                 nrow = 3, ncol = 4, byrow = TRUE))

index_dates <- list(matrix(c(1, 2), nrow = 2),
                    matrix(c(1, 2, 1,
                             2, 3, 4), nrow = 2, ncol = 3, byrow = TRUE))

test_that("compute_delta_group_delay_and_indiv handles delay in one group", {
  # In group 1, for individual 2, what is delay 1 (date 1 -> date 2)?
  delta <- compute_delta_group_delay_and_indiv(D = D, group_idx = 1,
                                               indiv_idx = 2, delay_idx = 1,
                                               index_dates = index_dates)
  # 25 - 12 = 13
  expect_equal(delta, 13)
})

test_that("compute_delta_group_delay_and_indiv selects correct delay in group 2", {
  # In group 2, for individual 3, what is delay 2 (date 2 -> date 3)?
  delta <- compute_delta_group_delay_and_indiv(D = D, group_idx = 2,
                                               indiv_idx = 3, delay_idx = 2,
                                               index_dates = index_dates)
  # 150 - 145 = 5
  expect_equal(delta, 5)
})

test_that("compute_delta_group_delay_and_indiv handles delay involving 4th date column", {
  # In group 2, for individual 1, what is delay 3 (date 1 -> date 4)?
  delta <- compute_delta_group_delay_and_indiv(D = D, group_idx = 2,
                                               indiv_idx = 1, delay_idx = 3,
                                               index_dates = index_dates)
  # 110 - 100 = 10
  expect_equal(delta, 10)
})

test_that("compute_delta_group_delay_and_indiv works for multiple individuals", {
  # In group 2, for individuals 2 and 3, what is delay 1 (date 1 -> date 2)?
  delta <- compute_delta_group_delay_and_indiv(D = D, group_idx = 2,
                                               indiv_idx = c(2, 3),
                                               delay_idx = 1,
                                               index_dates = index_dates)
  # indiv 2: 120 - 110 = 10; indiv 3: 145 - 130 = 15
  expect_equal(delta, c(10, 15))
})

# compute_delta --------------------------------------------------------------

D <- list(matrix(c(10, 20,
                   12, 25,
                   15, 28),
                 nrow = 3, ncol = 2, byrow = TRUE),
          matrix(c(100, 105, 108, 110,
                   110, 120, 125, 130,
                   130, 145, 150, 140),
                 nrow = 3, ncol = 4, byrow = TRUE))

index_dates <- list(matrix(c(1, 2), nrow = 2),
                    matrix(c(1, 2, 1,
                             2, 3, 4), nrow = 2, ncol = 3, byrow = TRUE))

test_that("compute_delta returns a list of the correct length", {
  delays <- compute_delta(D, index_dates)
  expect_true(is.list(delays))
  expect_length(delays, length(D))
})

test_that("compute_delta correctly calculates delay matrix for a simple group", {
  delays <- compute_delta(D, index_dates)
  # For group 1, index_dates has only one delay (date 1 -> date 2)
  # 20 - 10 = 10
  # 25 - 12 = 13
  # 28 - 15 = 13
  expected_g1 <- matrix(c(10, 13, 13))
  expect_equal(delays[[1]], expected_g1)
})

test_that("compute_delta correctly calculates delay matrix for multiple delays", {
  delays <- compute_delta(D, index_dates)
  # For group 2, index_dates has three delays
  # delay 1 (date 1 -> date 2): 105 - 100 = 5, 120 - 110 = 10, 145 - 130 = 15
  # delay 2 (date 2 -> date 3): 108 - 105 = 3, 125 - 120 = 5, 150 - 145 = 5
  # delay 3 (date 1 -> date 4): 110 - 100 = 10, 130 - 110 = 20, 140 - 130 = 10
  expected_matrix_g2 <- matrix(c(5, 3, 10,
                                 10, 5, 20,
                                 15, 5, 10),
                               nrow = 3, ncol = 3, byrow = TRUE)
  
  expect_equal(delays[[2]], expected_matrix_g2)
})

# find_range -----------------------------------------------------------------

test_that("find_range finds the correct min and max across 2 groups", {
  obs_dat <- list(matrix(c(10, 20, 15, 25), nrow = 2),
                  matrix(c(5, 50, 30, 40), nrow = 2))
  expected_range <- c(5, 50)
  expect_equal(find_range(obs_dat), expected_range)
})

test_that("find_range correctly handles NA values", {
  # NAs should be ignored when other dates available
  obs_dat_na <- list(matrix(c(10, NA, 15, 25), nrow = 2),
                     matrix(c(5, 50, NA, 40), nrow = 2))
  expect_equal(find_range(obs_dat_na), c(5, 50))
  
  # a date column in one group has all NAs
  obs_dat_na_col <- list(matrix(c(10, 20, NA, NA), nrow = 2),
                         matrix(c(5, 50, NA, 40), nrow = 2))
  expect_error(find_range(obs_dat_na_col),
    "All dates in column 2, group 1 are NA. Check that each date column in obs_dat has >=1 date."
  )
  
  # all values are NA
  obs_dat_all_na <- list(matrix(c(NA, NA, NA, NA), nrow = 2))
  expect_error(find_range(obs_dat_all_na),
    "All dates in column 1, group 1 are NA. Check that each date column in obs_dat has >=1 date."
  )
})

# compute_index_dates_order --------------------------------------------------

testthat::test_that("Simple transitivity is computed correctly", {
  # Output should include the transitive rule 1->3
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  out <- compute_index_dates_order(index_dates)
  expected <- list(cbind(c(1, 2), c(2, 3), c(1, 3)))
  testthat::expect_equal(out, expected)
})

testthat::test_that("Multiple transitive relationships computed correctly", {
  # Output should include 1->3, 2->4, and 1->4
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(3, 4)))
  out <- compute_index_dates_order(index_dates)
  expected <- list(cbind(c(1, 2), c(2, 3), c(3, 4), c(1, 3), c(2, 4), c(1, 4)))
  testthat::expect_equal(out, expected)
})

testthat::test_that("No transitive relationships handled correctly", {
  # No transitive relationships so output should be the same as input
  index_dates <- list(cbind(c(1, 2), c(3, 4)))
  out <- compute_index_dates_order(index_dates)
  testthat::expect_equal(out, index_dates)
})

testthat::test_that("Duplicates in the input are removed", {
  # Duplicate 1->2 should be removed (and transitive 1->3 should be added)
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(1, 2)))
  out <- compute_index_dates_order(index_dates)
  expected <- list(cbind(c(1, 2), c(2, 3), c(1, 3)))
  testthat::expect_equal(out, expected)
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

# check_MCMC_settings --------------------------------------------------------

index_dates <- list(
  matrix(c(1, 2), nrow = 2),
  cbind(c(1, 2), c(1, 3)),
  cbind(c(1, 2), c(2, 3),c(1, 4)),
  cbind(c(1, 2), c(2, 3), c(1, 4))
)

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
      0.25,
      c(0.25, 0.25),
      c(0.25, 0.25, 0.25),
      c(0.25, 0.25, 0.25))
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
