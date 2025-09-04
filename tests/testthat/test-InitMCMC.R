# Functions in InitMCMC.R
# - initialise_aug_data (x)
# - are_dates_incompatible (x)
# - initialise_theta_from_aug_data (x)

# initialise_aug_data --------------------------------------------------------

test_that("initialise_aug_data returns correctly structured output", {
  obs_dat <- list(matrix(c(10, 11), nrow = 1),
                  matrix(c(10, 15), nrow = 1))
  index_dates <- list(cbind(c(1, 2)),
                      cbind(c(1, 2)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 10))
  aug_dat <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  expect_type(aug_dat, "list")
  expect_named(aug_dat, c("D", "E"))
  expect_length(aug_dat$D, length(obs_dat))
  expect_equal(dim(aug_dat$D[[1]]), dim(obs_dat[[1]]))
  expect_equal(dim(aug_dat$D[[2]]), dim(obs_dat[[2]]))
})

test_that("Missing dates in obs_dat are imputed during initialisation and E is given -1", {
  obs_dat <- list(matrix(c(10, 16, NA,
                           12, 20, 27), nrow = 2, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 15))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # no NAs in output for D
  expect_false(any(is.na(out$D[[1]])))
  # missing date is given -1 in E
  expect_equal(out$E[[1]][1, 3], -1)
})

test_that("Dates involved in delays that are too long/short are modified", {
  # Make delays outside the bounds of mindelay = 3 and maxdelay = 10
  obs_dat <- list(matrix(c(10, 11, 20, # delay between 1->2 = 1 (< mindelay)
                           12, 20, 31 # delay between 2->3 = 11 (> maxdelay)
                           ), nrow = 2, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # expect either [1,1] or [1,2] and [2,2] or [2,3] to be flagged as an error
  expect_true(any(out$E[[1]][1, 1] == 1 | out$E[[1]][1, 2] == 1))
  expect_true(any(out$E[[1]][2, 2] == 1 | out$E[[1]][2, 3] == 1))
})

test_that("Date involved in >1 problematic delay is modified", {
  # date 2 for indiv 1 involved in 3 problematic delays
  obs_dat <- list(matrix(c(10, 80, 14, 18,
                           12, 20, 25, 27), nrow = 2, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(2, 4)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 15))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # expect date [1,2] to be modified and flagged as error
  expect_true(out$D[[1]][1, 2] != obs_dat[[1]][1, 2])
  expect_true(out$E[[1]][1, 2] == 1)
})

test_that("Date furthest from median is corrected when only one problematic delay exists", {
  # create only one problematic delay - delay 1->4 is 90 days (> maxdelay)
  obs_dat <- list(matrix(c(10, 15, 25, 100), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(1, 4)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 15))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # check that date furthest from the median (100) is corrected
  # and that the other date in the problematic delay (10) is not
  expect_true(out$D[[1]][1, 4] != obs_dat[[1]][1, 4])
  expect_true(out$D[[1]][1, 1] == obs_dat[[1]][1, 1])
  # check that [1, 4] is flagged as an error and the other date [1, 1] is not
  expect_equal(out$E[[1]][1, 4], 1)
  expect_equal(out$E[[1]][1, 1], 0)
})

test_that("Transitive error is corrected", {
  # transitive delay 1->4 is 70 days (> maxdelay * 3)
  obs_dat <- list(matrix(c(10, 15, 20, 80), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3), c(3, 4)))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 20))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # check that date furthest from the median (100) is corrected
  # and that the other date in the problematic delay (10) is not
  expect_true(out$D[[1]][1, 4] != obs_dat[[1]][1, 4])
  # check for correct error indicator
  expect_equal(out$E[[1]][1, 4], 1)
  expect_equal(out$E[[1]][1, 1], 0)
})

test_that("No dates are modified if there are no incompatible delays", {
  obs_dat <- list(matrix(c(10, 15, 20), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # output should be the same as the input
  expect_equal(out$D[[1]], obs_dat[[1]])
  # E should be all 0
  expect_equal(out$E[[1]][1, ], c(0, 0, 0))
})

test_that("Two groups are handled correctly", {
  # group 1 has an incompatible delay (1 day, mindelay = 3)
  # group 2 has a valid delay (5 days)
  obs_dat <- list(matrix(c(10, 11), nrow = 1),
                  matrix(c(10, 15), nrow = 1))
  index_dates <- list(cbind(c(1, 2)),
                      cbind(c(1, 2)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # check that exactly one error is flagged in group 1
  expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 1)
  # check that group 2 has no errors
  expect_true(all(out$E[[2]] == 0))
})

test_that("Mixed problematic groups are handled correctly", {
  # group 1 has a problematic delay (> maxdelay)
  # group 2 has a missing date
  obs_dat <- list(matrix(c(10, 50), nrow = 1),
                  matrix(c(10, NA), nrow = 1))
  index_dates <- list(cbind(c(1, 2)),
                      cbind(c(1, 2)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  # check that one error is flagged in group 1
  expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 1)
  # check that missing value in group 2 is imputed and flagged as -1
  expect_false(is.na(out$D[[2]][1, 2]))
  expect_equal(out$E[[2]][1, 2], -1)
})

test_that("All delays less than mindelay handled correcty", {
  # both specified delays are incompatible: 1->2 (delay = 1), 2->3 (delay = 2)
  # transitive delay also incompatible: 1->3 (delay = 2)
  obs_dat <- list(matrix(c(10, 11, 12), nrow = 1, byrow = TRUE))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 2)
})

# Edge case to discuss:
test_that("Specified delays invalid, transitive delay valid", {
  # both specified delays are incompatible: 1->2 (delay = 1), 2->3 (delay = 2)
  # under previous rules the transitive delay was compatible: 1->3 (delay = 3)
  # this has been changed (now 1->3 is mindelay * 2)
  obs_dat <- list(matrix(c(10, 11, 13), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 2)
  # Under previous rules date 2 is unchanged because although date 2 is
  # flagged as being in both bad pairs and set to NA, in the "missing" part
  # date 2 is constrained as both "before" date 3 and "after" date 1, so
  # floor(median(c(10, 13))) = 11, the same date as the original and so E = 0
})

test_that("All delays exceeding maxdelay handled correctly", {
  # both specified delays exceed maxdelay: 1->2 (delay = 11), 2->3 (delay = 11)
  # transitive delay also incompatible: 1->3 (delay = 22)
  obs_dat <- list(matrix(c(10, 21, 32), nrow = 1))
  index_dates <- list(cbind(c(1, 2), c(2, 3)))
  MCMC_settings <- list(init_options = list(mindelay = 3, maxdelay = 10))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  expect_equal(sum(out$E[[1]] == 1, na.rm = TRUE), 2)
})

test_that("Missing dates in multiple groups handled correctly", {
  index_dates <- list(matrix(c(1, 2), nrow = 2),
                      cbind(c(1, 2), c(1, 3)),
                      cbind(c(1, 2), c(2, 3), c(1, 4)),
                      cbind(c(1, 2), c(2, 3), c(1, 4)))
  obs_dat <- list(matrix(c(20, 24,
                           19, 23), nrow = 2, byrow = TRUE),
                  matrix(c(30, 35, NA,
                           NA, 23, 40), nrow = 2, byrow = TRUE),
                  matrix(c(17, 23, 27, NA,
                           27, 33, 41, 40), nrow = 2, byrow = TRUE),
                  matrix(c(23, 31, NA, 27,
                           8, 10, 28, 13), nrow = 2, byrow = TRUE))
  MCMC_settings <- list(init_options = list(mindelay = 1, maxdelay = 30))
  out <- initialise_aug_data(obs_dat, index_dates, MCMC_settings)
  
  expected_E <- list(matrix(0, nrow = 2, ncol = 2),
                     matrix(c(0, 0, -1, -1, 0, 0), nrow = 2, byrow = TRUE),
                     matrix(c(0, 0, 0, -1, 0, 0, 0, 0), nrow = 2, byrow = TRUE),
                     matrix(c(0, 0, -1, 0, 0, 0, 0, 0), nrow = 2, byrow = TRUE))
  expect_equal(out$E, expected_E)
})

# are_dates_incompatible -----------------------------------------------------

test_that("are_dates_incompatible works correctly", {
  mindelay <- 3
  maxdelay <- 10
  expect_true(are_dates_incompatible(10, 11, mindelay, maxdelay)) # too short
  expect_true(are_dates_incompatible(11, 13, mindelay, maxdelay)) # too short
  expect_false(are_dates_incompatible(10, 13, mindelay, maxdelay)) # mindelay
  expect_false(are_dates_incompatible(10, 20, mindelay, maxdelay)) # maxdelay
  expect_true(are_dates_incompatible(10, 21, mindelay, maxdelay)) # too long
})

# initialise_theta_from_aug_data ---------------------------------------------

test_that("initialise_theta_from_aug_dat returns correctly structured output", {
  aug_dat <- list(
    D = list(matrix(c(10, 20,
                      12, 25), nrow = 2, ncol = 2, byrow = TRUE),
             matrix(c(100, 110, 130,
                      105, 120, 145), nrow = 2, ncol = 3, byrow = TRUE))
    )
  index_dates <- list(matrix(c(1, 2), nrow = 2),
                      matrix(c(1, 2, 2, 3), nrow = 2, byrow = TRUE))
  theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)
  
  expect_type(theta, "list")
  expect_named(theta, c("zeta", "mu", "CV"))
  expect_length(theta$zeta, 1)
  expect_length(theta$mu, length(aug_dat$D))
  expect_length(theta$CV, length(aug_dat$D))
  expect_length(theta$mu[[1]], ncol(index_dates[[1]]))
  expect_length(theta$mu[[2]], ncol(index_dates[[2]]))
  expect_length(theta$CV[[1]], ncol(index_dates[[1]]))
  expect_length(theta$CV[[2]], ncol(index_dates[[2]]))
})

test_that("zeta is initialised to the correct value", {
  aug_dat <- list(D = list(matrix(c(10, 20,
                                    15, 21), nrow = 2, byrow = TRUE)))
  index_dates <- list(matrix(c(1, 2), nrow = 2))
  
  zeta_default <- initialise_theta_from_aug_dat(aug_dat, index_dates)
  expect_equal(zeta_default$zeta, 0.1)

  zeta_custom <- initialise_theta_from_aug_dat(aug_dat, index_dates,
                                               zeta_init = 0.5)
  expect_equal(zeta_custom$zeta, 0.5)
})

test_that("mu and CV are calculated correctly for multiple groups", {
  aug_dat <- list(D = list(matrix(c(10, 20,
                                    12, 25,
                                    15, 28),
                                  nrow = 3, ncol = 2, byrow = TRUE),
                           matrix(c(100, 110, 130,
                                    105, 120, 145),
                                  nrow = 2, ncol = 3, byrow = TRUE)))
  index_dates <- list(matrix(c(1, 2), nrow = 2),
                      matrix(c(1, 2, 2, 3), nrow = 2, byrow = TRUE))
  theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)
  
  # group 1 delays: 10->20, 12->25, 15->28 -> c(10, 13, 13)
  delays_g1 <- c(10, 13, 13)
  expected_mu1 <- mean(delays_g1)
  expected_cv1 <- sd(delays_g1) / expected_mu1
  expect_equal(theta$mu[[1]], expected_mu1)
  expect_equal(theta$CV[[1]], expected_cv1)
  
  # group 2 delay 1: 100->110, 105->120 -> c(10, 15)
  # group 2 delay 2: 110->130, 120->145 -> c(20, 25)
  delays_g2_d1 <- c(10, 15)
  delays_g2_d2 <- c(20, 25)
  expected_mu2 <- c(mean(delays_g2_d1), mean(delays_g2_d2))
  expected_cv2 <- c(sd(delays_g2_d1) / mean(delays_g2_d1), 
                    sd(delays_g2_d2) / mean(delays_g2_d2))
  expect_equal(theta$mu[[2]], expected_mu2)
  expect_equal(theta$CV[[2]], expected_cv2)
})

test_that("initialise_theta_from_aug_dat handles groups with 1 individual", {
  aug_dat <- list(D = list(matrix(c(10, 25), nrow = 1)))
  index_dates <- list(matrix(c(1, 2), nrow = 2))
  theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)
  
  # mean of single delay (15) is just 15
  expect_equal(theta$mu[[1]], 15)
  # sd() of a single number is NA, so CV should be NA
  expect_true(is.na(theta$CV[[1]]))
})
