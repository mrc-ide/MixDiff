# Functions in LikelihoodPrior.R
# - lprior_total (x)
# - lprior_prob_error (x)
# - lprior_params_delay (x)
# - compute_n_errors (x)
# - lposterior_total (x)
# - LL_total (x)
# - LL_observation_term (x)
# - LL_observation_term_by_group_delay_and_indiv (x)
# - LL_error_term (x)
# - LL_error_term_by_group_delay_and_indiv (x)
# - LL_delays_term (x)
# - LL_delays_term_by_group_delay_and_indiv (x)

#------------------------------------------------------------------------------
# Priors
#------------------------------------------------------------------------------

# lprior_prob_error() ---------------------------------------------------------

test_that("lprior_prob_error correcty calculates the log-prior for zeta", {
  theta <- list(zeta = 0.1)
  hyperparameters <- list(shape1_prob_error = 1, shape2_prob_error = 10)
  out <- lprior_prob_error(theta, hyperparameters)
  exp_out <- dbeta(theta$zeta,
                   hyperparameters$shape1_prob_error,
                   hyperparameters$shape2_prob_error,
                   log = TRUE)
  expect_equal(out, exp_out)
})

test_that("lprior_prob_error returns -Inf for zeta outside [0, 1]", {
  hyperparameters <- list(shape1_prob_error = 1, shape2_prob_error = 10)
  
  # zeta > 1 should be impossible
  theta_high <- list(zeta = 1.1)
  expect_equal(lprior_prob_error(theta_high, hyperparameters), -Inf)
  
  # zeta < 0 should be impossible
  theta_low <- list(zeta = -0.1)
  expect_equal(lprior_prob_error(theta_low, hyperparameters), -Inf)
})

# lprior_params_delay() -------------------------------------------------------

test_that("lprior_params_delay correctly calculates log-prior for mu and cv", {
  theta <- list(zeta = 0.1,
                mu = list(12, c(5, 10)),
                CV = list(0.15, c(0.3, 0.2)))
  hyperparameters <- list(shape1_prob_error = 1, shape2_prob_error = 10,
                          mean_mean_delay = 15, mean_CV_delay = 1)
  all_mu <- unlist(theta$mu)
  mu_out <- lprior_params_delay("mu", theta, hyperparameters)
  mu_exp <- sum(dexp(all_mu,
                     rate = 1 / hyperparameters$mean_mean_delay,
                     log = TRUE))
  expect_equal(mu_out, mu_exp)
  
  all_cv <- unlist(theta$CV)
  cv_out <- lprior_params_delay("CV", theta, hyperparameters)
  cv_exp <- sum(dexp(all_cv,
                     rate = 1 / hyperparameters$mean_mean_delay,
                     log = TRUE))
  expect_equal(cv_out, cv_exp)
})

# lprior_total() -------------------------------------------------------

test_that("lprior_total correctly sums individual log priors", {
  theta <- list(zeta = 0.1,
                mu = list(12, c(5, 10)),
                CV = list(0.15, c(0.3, 0.2)))
  hyperparameters <- list(shape1_prob_error = 1, shape2_prob_error = 10,
                          mean_mean_delay = 15, mean_CV_delay = 1)
  prior_zeta <- lprior_prob_error(theta, hyperparameters)
  prior_mu <- lprior_params_delay("mu", theta, hyperparameters)
  prior_cv <- lprior_params_delay("CV", theta, hyperparameters)
  
  out <- lprior_total(theta, hyperparameters)
  exp <- prior_zeta + prior_mu + prior_cv
  
  expect_equal(out, exp)
})

#------------------------------------------------------------------------------
# Observation terms
#------------------------------------------------------------------------------

# LL_observation_term_by_group_delay_and_indiv() ------------------------------

test_that("LL_observation_term_by_group_delay_and_indiv works correctly", {
  obs_dat <- list(matrix(c(10, 15, 26,
                           11, NA, 30),
                         nrow = 2, byrow = TRUE))
  
  aug_dat <- list(D = list(matrix(c(10, 15, 25, # date 3 is actually 25
                                    11, 19, 30),  # date 2 is actually 19
                                  nrow = 2, byrow = TRUE)),
                  E = list(matrix(c(0, 0, 1, # date 3 recorded with error
                                    0, -1, 0), # date 2 was missing
                                  nrow = 2, byrow = TRUE)))
  
  range_dates <- c(1, 50)
  
  # indiv 1, date 1 = no error, dates in aug_dat and obs_dat match -> log(1) = 0
  out_match <- LL_observation_term_by_group_delay_and_indiv(
    aug_dat, obs_dat, group_idx = 1, date_idx = 1, indiv_idx = 1, range_dates)
  expect_equal(out_match[1,1], 0)
  
  # indiv 1, date 3 = error (dates don't match) -> log(1 / (50-1))
  # assumes true date could be any date within range_dates with equal prob
  # K = 1 / diff(range_dates)
  out_error <- LL_observation_term_by_group_delay_and_indiv(
    aug_dat, obs_dat, group_idx = 1, date_idx = 3, indiv_idx = 1, range_dates)
  expect_equal(out_error[1,1], log(1 / as.numeric(diff(range_dates))))
  
  # indiv 2, date 2 = missing -> log(1 / (50-1))
  # same as above - date could be any date within range_dates with equal prob
  out_missing <- LL_observation_term_by_group_delay_and_indiv(
    aug_dat, obs_dat, group_idx = 1, date_idx = 2, indiv_idx = 2, range_dates)
  expect_equal(out_missing[1,1], log(1 / as.numeric(diff(range_dates))))
  
  # sanity check - no error, but data mismatch -> log(0) = -1e5
  # this shouldn't be possible, if there's no error then aug_dat = obs_dat
  # to avoid problems with -Inf, it is replaced by a large negative number
  aug_dat_mismatch <- aug_dat
  aug_dat_mismatch$D[[1]][1, 1] <- 99 # force a mismatch
  out_mismatch <- LL_observation_term_by_group_delay_and_indiv(
    aug_dat_mismatch, obs_dat, group_idx = 1, date_idx = 1, indiv_idx = 1,
    range_dates)
  expect_equal(out_mismatch[1, 1], -1e5)
})

test_that("LL_observation_term correctly sums individual likelihoods", {
  obs_dat <- list(matrix(c(10, 15, 26,
                           11, NA, 30),
                         nrow = 2, byrow = TRUE))
  aug_dat <- list(D = list(matrix(c(10, 15, 25, # date 3 is actually 25
                                    11, 19, 30),  # date 2 is actually 19
                                  nrow = 2, byrow = TRUE)),
                  E = list(matrix(c(0, 0, 1, # date 3 recorded with error
                                    0, -1, 0), # date 2 was missing
                                  nrow = 2, byrow = TRUE)))
  range_dates <- c(1, 50)
  out <- LL_observation_term(aug_dat, obs_dat, range_dates)
  # sum of all individual observation likelihoods
  manual_sum <- sum(LL_observation_term_by_group_delay_and_indiv(
    aug_dat, obs_dat, group_idx = 1, 
    date_idx = seq_len(ncol(aug_dat$D[[1]])), 
    indiv_idx = seq_len(nrow(aug_dat$D[[1]])), 
    range_dates
  ))
  expect_equal(out, manual_sum)
})


#------------------------------------------------------------------------------
# Error terms
#------------------------------------------------------------------------------

# compute_n_errors() ----------------------------------------------------------

test_that("compute_n_errors correctly counts errors and recorded dates", {
  aug_dat <- list(E = list(matrix(c( 1,  0, -1,
                                     0,  1,  1), nrow = 2, byrow = TRUE)))
  result <- compute_n_errors(aug_dat)
  expect_equal(result[1], 3) # number of errors
  expect_equal(result[2], 5) # number of recorded dates
})

# LL_error_term_by_group_delay_and_indiv() ------------------------------------

test_that("LL_error_term_by_group_delay_and_indiv works correctly", {
  obs_dat <- list(matrix(c(10, 15, 26,
                           11, NA, 30),
                         nrow = 2, byrow = TRUE))
  aug_dat <- list(D = list(matrix(c(10, 15, 25, # date 3 is actually 25
                                    11, 19, 30),  # date 2 is actually 19
                                  nrow = 2, byrow = TRUE)),
                  E = list(matrix(c(0, 0, 1, # date 3 recorded with error
                                    0, -1, 0), # date 2 was missing
                                  nrow = 2, byrow = TRUE)))
  theta <- list(
    zeta = 0.1, # 10% chance of error
    mu = list(c(5, 10)), # Mean delays
    CV = list(c(0.2, 0.3))  # CV for delays
  )
  
  # indiv 1, date 1 is E = 0 (no error) -> log(1 - zeta)
  out_correct <- LL_error_term_by_group_delay_and_indiv(
    aug_dat, theta, group_idx = 1, date_idx = 1, indiv_idx = 1)
  expect_equal(out_correct[1, 1], log(1 - theta$zeta))
  
  # indiv 1, date 3 is E = 1 (error) -> log(zeta)
  out_error <- LL_error_term_by_group_delay_and_indiv(
    aug_dat, theta, group_idx = 1, date_idx = 3, indiv_idx = 1)
  expect_equal(out_error[1, 1], log(theta$zeta))
  
  # indiv 2, date 2 is E = -1 (missing) -> 0
  out_missing <- LL_error_term_by_group_delay_and_indiv(
    aug_dat, theta, group_idx = 1, date_idx = 2, indiv_idx = 2)
  expect_equal(out_missing[1, 1], 0)
})

# LL_error_term() -------------------------------------------------------------

test_that("LL_error_term returns expected log-prob with known errors", {
  aug_dat <- list(E = list(matrix(c(1, 0, -1, 1), nrow = 2)))
  theta <- list(zeta = 0.25)
  result <- LL_error_term(aug_dat, theta)
  # 2 errors, 1 correct, 1 missing -> 3 recorded values in total
  # 2 recorded with error = log(0.25) * 2, 1 recorded correctly = log(0.75) * 1
  expected <- log(0.25) * 2 + log(0.75) * 1
  expect_equal(result, expected)
})

test_that("LL_error_term aggregates errors across multiple groups", {
  aug_dat <- list(E = list(matrix(c(1, 0, 0, 0), nrow = 2), 
                           matrix(c(1, 1, 0, 0), nrow = 2)))
  theta <- list(zeta = 0.25)
  # 3 errors, 5 correct
  result <- LL_error_term(aug_dat, theta)
  expected <- log(0.25) * 3 + log(1 - 0.25) * 5
  expect_equal(result, expected)
})


#------------------------------------------------------------------------------
# Delay terms
#------------------------------------------------------------------------------

# LL_delays_term_by_group_delay_and_indiv() -----------------------------------

test_that("LL_delays_term_by_group_delay_and_indiv handles one individual", {
  aug_dat <- list(D = list(matrix(c(10, 20,
                                    12, 25), nrow = 2, byrow = TRUE)))
  obs_dat <- aug_dat$D
  index_dates <- list(matrix(c(1, 2), nrow = 2))
  theta <- list(mu = list(10), CV = list(0.15))
  
  ll_out <- LL_delays_term_by_group_delay_and_indiv(
    aug_dat, theta, obs_dat, 
    group_idx = 1, delay_idx = 1, indiv_idx = 2, 
    index_dates = index_dates
  )

  # for group 1, indiv 1, the delay is 25 - 12 = 13
  expected_ll <- DiscrGamma(k = 13, 
                            mu = theta$mu[[1]], 
                            cv = theta$CV[[1]], 
                            log = TRUE)
  
  expect_equal(ll_out, expected_ll)
})

test_that("LL_delays_term_by_group_delay_and_indiv handles multiple individuals", {
  aug_dat <- list(D = list(matrix(c(100, 110, 130,
                                    105, 120, 145), nrow = 2, byrow = TRUE)))
  obs_dat <- aug_dat$D
  index_dates <- list(matrix(c(2, 3), nrow = 2))
  theta <- list(mu = list(20), CV = list(0.15))
  
  ll_out <- LL_delays_term_by_group_delay_and_indiv(
    aug_dat, theta, obs_dat, 
    group_idx = 1, delay_idx = 1, indiv_idx = c(1, 2), 
    index_dates = index_dates
  )
  
  # delays are: 130 - 110 = 20 (indiv 1), 145 - 120 = 25 (indiv 2)
  delays <- c(20, 25)
  expected_ll_vec <- DiscrGamma(k = delays, 
                                mu = theta$mu[[1]], 
                                cv = theta$CV[[1]], 
                                log = TRUE)
  
  expect_equal(ll_out, expected_ll_vec)
})

# LL_delays_term() ------------------------------------------------------------

test_that("LL_delays_term correctly sums individual likelihoods", {
  aug_dat <- list(D = list(matrix(c(10, 20,
                                    12, 25), nrow = 2, byrow = TRUE),
                           matrix(c(100, 110, 130, 135,
                                    105, 120, 145, 140), nrow = 2, byrow = TRUE)))
  obs_dat <- aug_dat$D
  index_dates <- list(matrix(c(1, 2), nrow = 2),
                      matrix(c(1, 2, 1,
                               2, 3, 4), nrow = 2, byrow = TRUE))
  theta <- list(mu = list(11.5, c(12.5, 22.5, 40)),
                CV = list(0.15, c(0.15, 0.15, 0.15)))
  
  total_ll_out <- LL_delays_term(aug_dat, theta, obs_dat, index_dates)
  
  # manually calculate the total by summing the parts
  ll_g1_d1 <- sum(DiscrGamma(c(10, 13), theta$mu[[1]], theta$CV[[1]], log = TRUE))
  ll_g2_d1 <- sum(DiscrGamma(c(10, 15), theta$mu[[2]][1], theta$CV[[2]][1], log = TRUE))
  ll_g2_d2 <- sum(DiscrGamma(c(20, 25), theta$mu[[2]][2], theta$CV[[2]][2], log = TRUE))
  ll_g2_d3 <- sum(DiscrGamma(c(35, 35), theta$mu[[2]][3], theta$CV[[2]][3], log = TRUE))
  expected_total_ll <- ll_g1_d1 + ll_g2_d1 + ll_g2_d2 + ll_g2_d3
  
  expect_equal(total_ll_out, expected_total_ll)
})

test_that("LL_delays_term returns a single numeric value", {
  aug_dat <- list(D = list(matrix(c(10, 20, 12, 25), nrow = 2, byrow = TRUE)))
  obs_dat <- aug_dat$D
  index_dates <- list(matrix(c(1, 2), nrow = 2))
  theta <- list(mu = list(11), CV = list(0.15))
  
  ll_out <- LL_delays_term(aug_dat, theta, obs_dat, index_dates)
  
  expect_true(is.numeric(ll_out))
  expect_length(ll_out, 1)
})

#------------------------------------------------------------------------------
# Total
#------------------------------------------------------------------------------

# LL_total() ------------------------------------------------------------------

test_that("LL_total correctly sums all likelihood components", {
  aug_dat <- list(D = list(matrix(c(10, 20,
                                    12, 25), nrow = 2, byrow = TRUE),
                           matrix(c(100, 110, 130, 135,
                                    105, 120, 145, 140), nrow = 2, byrow = TRUE)),
                  E = list(matrix(c(0, 1,
                                    0, 0), nrow = 2, byrow = TRUE),
                           matrix(c(0, 1, -1, 0,
                                    0, 0, 0, 1), nrow = 2, byrow = TRUE)))
  obs_dat <- list(matrix(c(10, 50,
                           12, 25), nrow = 2, byrow = TRUE),
                  matrix(c(100, 80, NA, 135,
                           105, 120, 145, 200), nrow = 2, byrow = TRUE))
  
  index_dates <- list(matrix(c(1, 2), nrow = 2),
                      matrix(c(1, 2, 1,
                               2, 3, 4), nrow = 2, byrow = TRUE))
  theta <- list(zeta = 0.1,
                mu = list(11.5, c(12.5, 22.5, 40)),
                CV = list(0.15, c(0.15, 0.15, 0.15)))
  range_dates <- c(1, 150)
  
  out <- LL_total(aug_dat, theta, obs_dat, index_dates, range_dates)
  
  exp_out <- LL_observation_term(aug_dat, obs_dat, range_dates) +
    LL_error_term(aug_dat, theta) +
    LL_delays_term(aug_dat, theta, obs_dat, index_dates)
  
  expect_equal(out, exp_out)
})


# lposterior_total() -------------------------------------------------------

test_that("lposterior_total correctly sums likelihood and prior", {
  aug_dat <- list(D = list(matrix(c(10, 20,
                                    12, 25), nrow = 2, byrow = TRUE),
                           matrix(c(100, 110, 130, 135,
                                    105, 120, 145, 140), nrow = 2, byrow = TRUE)),
                  E = list(matrix(c(0, 1,
                                    0, 0), nrow = 2, byrow = TRUE),
                           matrix(c(0, 1, -1, 0,
                                    0, 0, 0, 1), nrow = 2, byrow = TRUE)))
  obs_dat <- list(matrix(c(10, 50,
                           12, 25), nrow = 2, byrow = TRUE),
                  matrix(c(100, 80, NA, 135,
                           105, 120, 145, 200), nrow = 2, byrow = TRUE))
  
  index_dates <- list(matrix(c(1, 2), nrow = 2),
                      matrix(c(1, 2, 1,
                               2, 3, 4), nrow = 2, byrow = TRUE))
  theta <- list(zeta = 0.1,
                mu = list(11.5, c(12.5, 22.5, 40)),
                CV = list(0.15, c(0.15, 0.15, 0.15)))
  range_dates <- c(1, 150)
  hyperparameters <- list(shape1_prob_error = 1, shape2_prob_error = 10,
                          mean_mean_delay = 15, mean_CV_delay = 1)
  
  out <- lposterior_total(aug_dat, theta, obs_dat,
                          hyperparameters, index_dates, range_dates)
  
  exp_out <- LL_total(aug_dat, theta, obs_dat, index_dates, range_dates) +
    lprior_total(theta, hyperparameters)
  
  expect_equal(out, exp_out)
})
