# LL_observation_term_by_group_delay_and_indiv() ------------------------------

test_that("LL_observation_term_by_group_delay_and_indiv handles E = 0", {
  aug_dat <- list(
    D = list(matrix(c(1, 2, 3), nrow = 3, ncol = 1)),
    E = list(matrix(0, nrow = 3, ncol = 1))
  )
  obs_dat <- list(matrix(c(1, 2, 3), nrow = 3, ncol = 1))
  theta <- list(zeta = 0.1)
  
  ll <- LL_observation_term_by_group_delay_and_indiv(
    aug_dat, theta, obs_dat, group_idx = 1,
    date_idx = 1, indiv_idx = 1:3,
    range_dates = c(1, 10)
  )
  
  expect_true(all(ll == 0))
})

test_that("LL_observation_term_by_group_delay_and_indiv handles E = 1", {
  aug_dat <- list(
    D = list(matrix(c(5, 5, 5), nrow = 3, ncol = 1)),
    E = list(matrix(1, nrow = 3, ncol = 1))
  )
  obs_dat <- list(matrix(c(9, 9, 9), nrow = 3, ncol = 1))
  theta <- list(zeta = 0.1)
  
  ll <- LL_observation_term_by_group_delay_and_indiv(
    aug_dat, theta, obs_dat, group_idx = 1,
    date_idx = 1, indiv_idx = 1:3,
    range_dates = c(1, 10)
  )
  
  expect_true(all(is.finite(ll)))
  # If the true date is incorrectly observed then any date in the valid range
  # is equally likely (uniform probability of drawing any date in range_dates)
  # K = 1 / diff(range_dates) = 1/9
  expect_true(all(ll == log(1/9)))
})

# compute_n_errors() ----------------------------------------------------------

test_that("compute_n_errors correctly counts errors and recorded dates", {
  aug_dat <- list(E = list(
    matrix(c( 1,  0, -1,
              0,  1,  1), nrow = 2)
  ))
  
  obs_dat <- NULL
  
  result <- compute_n_errors(aug_dat, obs_dat)
  
  expect_equal(result[1], 3) # number of errors
  expect_equal(result[2], 5) # number of recorded dates
})

# LL_error_term() -------------------------------------------------------------

test_that("LL_error_term returns expected log-prob with known errors", {
  aug_dat <- list( E = list(matrix(c(1, 0, -1, 1), nrow = 2)))
  obs_dat <- list(matrix(1:4, nrow = 2))
  theta <- list(zeta = 0.25)
  
  result <- LL_error_term(aug_dat, theta, obs_dat)
  
  expected <- log(0.25) * 2 + log(0.75) * 1  # 2 errors, 1 correct, 1 missing
  expect_equal(result, expected)
})


# LL_delays_term_by_group_delay_and_indiv() -----------------------------------

test_that("LL_delays_term_by_group_delay_and_indiv returns valid log-density", {
  aug_dat <- list(
    D = list(matrix(c(1, 6), ncol = 2))
  )
  obs_dat <- list(matrix(NA, nrow = 1, ncol = 2))
  theta <- list(mu = list(5), CV = list(0.5))
  index_dates <- list(matrix(c(1, 2), nrow = 2))
  
  log_ll <- LL_delays_term_by_group_delay_and_indiv(
    aug_dat, theta, obs_dat, group_idx = 1,
    delay_idx = 1, indiv_idx = 1,
    index_dates = index_dates
  )
  
  expect_type(log_ll, "double")
  expect_true(is.finite(log_ll))
})

# LL_total() ------------------------------------------------------------------

test_that("LL_total returns finite value for simple valid input", {
  theta <- list(mu = list(5), CV = list(0.5), zeta = 0.1)
  index_dates <- list(matrix(c(1, 2), nrow = 2))
  
  aug_dat <- list(
    D = list(matrix(c(1, 6), ncol = 2, byrow = TRUE)),
    E = list(matrix(c(0, 1), ncol = 2))
  )
  obs_dat <- list(matrix(c(1, 3), ncol = 2))
  
  result <- LL_total(aug_dat, theta, obs_dat, index_dates, range_dates = c(1, 10))
  
  expect_type(result, "double")
  expect_true(is.finite(result))
})

# lprior_prob_error() ---------------------------------------------------------

test_that("find_params_beta returns valid beta parameters", {
  param_beta <- find_params_beta(mean = 0.2, var = 0.01)

  # Check structure
  expect_type(param_beta, "double")
  expect_length(param_beta, 2)
  expect_true(all(param_beta > 0))
  
})

test_that("lprior_prob_error works as expected", {
  param_beta <- find_params_beta(mean = 0.2, var = 0.01)
  theta <- list(zeta = 0.2)
  hyperparams <- list(shape1_prob_error = param_beta[1],
                      shape2_prob_error = param_beta[2])
  
  log_prior <- lprior_prob_error(theta, hyperparams)
  
  expect_type(log_prior, "double")
  expect_true(is.finite(log_prior))
})

# lprior_params_delay() -------------------------------------------------------

test_that("lprior_params_delay returns finite log-density for valid delays", {
  theta <- list(
    mu = list(c(5, 10), c(7, 8)),
    CV = list(c(0.5, 0.6), c(0.4, 0.7))
  )
  hyperparams <- list(mean_mean_delay = 10, mean_CV_delay = 10)
  
  log_prior_mu <- lprior_params_delay("mu", theta, hyperparams)
  log_prior_cv <- lprior_params_delay("CV", theta, hyperparams)
  
  expect_type(log_prior_mu, "double")
  expect_type(log_prior_cv, "double")
  expect_true(is.finite(log_prior_mu))
  expect_true(is.finite(log_prior_cv))
})

# lprior_total() -------------------------------------------------------

test_that("lprior_total combines priors correctly and returns finite value", {
  theta <- list(
    mu = list(c(5, 10), c(7, 8)),
    CV = list(c(0.5, 0.6), c(0.4, 0.7)),
    zeta = 0.1
  )
  
  hyperparams <- list(
    shape1_prob_error = 3,
    shape2_prob_error = 12,
    mean_mean_delay = 100,
    mean_CV_delay = 100
  )
  
  log_prior <- lprior_total(theta, hyperparams)
  
  expect_type(log_prior, "double")
  expect_true(is.finite(log_prior))
})

# lposterior_total() -------------------------------------------------------

test_that("lposterior_total returns finite value for simulated data", {

  set.seed(10)
  theta <- list(
    prop_missing_data = 0.2,
    zeta = 0.05,
    mu = list(5, c(10, 15)),
    CV = list(0.5, c(0.5, 0.5))
  )
  n_groups <- 2
  n_per_group <- c(10, 10)
  range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
                               as.Date("01/01/2015", "%d/%m/%Y")))
  index_dates <- list(matrix(c(1, 2), nrow = 2),
                      cbind(c(1, 2), c(1, 3)))
  
  D <- simul_true_data(theta, n_per_group, range_dates, index_dates,
                       simul_error = TRUE)
  
  aug_dat <- list(D = D$true_dat, E = D$E)
  obs_dat <- D$obs_dat
  
  hyperparams <- list(
    shape1_prob_error = 3,
    shape2_prob_error = 12,
    mean_mean_delay = 100,
    mean_CV_delay = 100
  )
  
  log_post <- lposterior_total(aug_dat, theta, obs_dat, hyperparams, index_dates)
  
  expect_type(log_post, "double")
  expect_true(is.finite(log_post))
})
