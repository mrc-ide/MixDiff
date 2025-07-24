# Test simul_true_data() and simul_obs_data()

set.seed(1)

# Create data
n_groups <- 4
n_dates <- c(2, 3, 4, 4)

mu_baseline <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
cv_baseline <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))

theta <- list(
  prop_missing_data = 0.2,
  zeta = 0.05,
  mu = mu_baseline,
  CV = cv_baseline
)

n_per_group <- rep(100, n_groups)

range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
                             as.Date("01/01/2015", "%d/%m/%Y")))

index_dates <- list(
  matrix(c(1, 2), nrow = 2),
  cbind(c(1, 2), c(1, 3)),
  cbind(c(1, 2), c(2, 3),c(1, 4)),
  cbind(c(1, 2), c(2, 3), c(1, 4))
)

index_dates_order <- list(
  matrix(c(1, 2), nrow = 2),
  cbind(c(1, 2), c(1, 3)),
  cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)),
  cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4))
)

# ----------------------------------------------------------------------------

test_that("structure of simul_true_data output is correct", {
  true_dataset <- simul_true_data(theta,
                                  n_per_group,
                                  range_dates,
                                  index_dates,
                                  simul_error = TRUE,
                                  remove_allNA_indiv = TRUE)
  
  # Test list structure
  expect_type(true_dataset, "list")
  expect_named(true_dataset, c("true_dat", "obs_dat", "E"))
  expect_length(true_dataset$true_dat, n_groups)
  expect_length(true_dataset$obs_dat, n_groups)
  expect_length(true_dataset$E, n_groups)
  
  # Test matrix structure
  # - nrow should be number of individuals in each group
  # - ncol should be number of dates in each group
  for (g in seq_len(n_groups)) {
    expect_equal(nrow(true_dataset$true_dat[[g]]), n_per_group[g],
                 label = paste("true_dat group", g, "row count mismatch"))
    expect_equal(ncol(true_dataset$true_dat[[g]]), n_dates[g],
                 label = paste("true_dat group", g, "col count mismatch"))
  }
})

test_that("dimensions of true_dat, obs_dat and E are the same", {
  true_dataset <- simul_true_data(theta,
                                  n_per_group,
                                  range_dates,
                                  index_dates,
                                  simul_error = TRUE,
                                  remove_allNA_indiv = TRUE)
  
  for (g in seq_len(n_groups)) {
    expect_equal(dim(true_dataset$true_dat[[g]]), dim(true_dataset$obs_dat[[g]]))
    expect_equal(dim(true_dataset$obs_dat[[g]]), dim(true_dataset$E[[g]]))
  }
})

# ----------------------------------------------------------------------------

test_that("obs_dat and E are NULL is simul_error is FALSE", {
  
  true_dataset <- simul_true_data(
    theta,
    n_per_group,
    range_dates,
    index_dates,
    simul_error = FALSE
  )
  
  expect_null(true_dataset$obs_dat)
  expect_null(true_dataset$E)
  
})

# ----------------------------------------------------------------------------

test_that("mean delays in true_dat match theta$mu", {
  set.seed(1)
  
  true_dataset <- simul_true_data(
    theta,
    n_per_group,
    range_dates,
    index_dates,
    simul_error = FALSE
  )
  
  tolerance <- 0.1  # 0.1 days
  
  for (g in seq_along(true_dataset$true_dat)) {
    dat <- true_dataset$true_dat[[g]]
    idx <- index_dates[[g]]
    mu_vec <- theta$mu[[g]]
    # dat <- true_dataset$true_dat[[3]]
    # idx <- index_dates[[3]]
    # mu_vec <- theta$mu[[3]]
    
    for (j in seq_len(ncol(idx))) {
      origin <- idx[1, j]
      target <- idx[2, j]
      # origin <- idx[1, 1]
      # target <- idx[2, 1]
      # origin <- idx[1, 2]
      # target <- idx[2, 2]
      # origin <- idx[1, 3]
      # target <- idx[2, 3]

      # Compute observed delays
      delay <- dat[, target] - dat[, origin]
#browser()
      # Check mean delay approx theta$mu[[g]][j]
      expect_true(
        abs(mean(delay) - mu_vec[j]) < tolerance,
        label = paste("Group", g, "Delay", j,
                      "mean =", round(mean(delay), 2),
                      "expected =", mu_vec[j])
      )
    }
  }
})

# ----------------------------------------------------------------------------

test_that("discr_gamma_sample works as expected", {

  delay_sim <- discr_gamma_sample(100, mu = 5, cv = 0.5)
  expect_equal(mean(delay_sim), 5, tolerance = 0.1)
  
  delay_sim <- discr_gamma_sample(100, mu = 7, cv = 0.5)
  expect_equal(mean(delay_sim), 7, tolerance = 0.1)
  
  delay_sim <- discr_gamma_sample(100, mu = 10, cv = 0.5)
  expect_equal(mean(delay_sim), 10, tolerance = 0.1)
  
  delay_sim <- discr_gamma_sample(100, mu = 12, cv = 0.5)
  expect_equal(mean(delay_sim), 12, tolerance = 0.1)

})

# ----------------------------------------------------------------------------

test_that(
  "correct prop of dates are missing in obs_dat when n_per_group = 100", {
  # smaller group
  n_per_group <- rep(100, n_groups)
  
  true_dataset <- simul_true_data(
    theta,
    n_per_group,
    range_dates,
    index_dates,
    simul_error = TRUE
  )
  
  target_missing_prop <- theta$prop_missing_data # 0.2
  tolerance <- 0.1  # 0.1 = between 10-30%
  
  for (g in seq_along(true_dataset$obs_dat)) {
    dat <- true_dataset$obs_dat[[g]]
    for (j in seq_len(ncol(dat))) {
      na_prop <- mean(is.na(dat[, j]))

      # Expect greater than 0.2 - tolerance
      expect_gt(
        na_prop,
        target_missing_prop - tolerance,
        label = paste("Group", g, "column", j, "has too few NAs:", na_prop)
      )
      
      # Expect less than 0.2 + tolerance
      expect_lt(
        na_prop,
        target_missing_prop + tolerance,
        label = paste("Group", g, "column", j, "has too many NAs:", na_prop)
      )
    }
  }
})
  
test_that(
  "correct prop of dates are missing in obs_dat when n_per_group = 1000", {   
  # larger group
  n_per_group <- rep(1000, n_groups)
  
  true_dataset <- simul_true_data(
    theta,
    n_per_group,
    range_dates,
    index_dates,
    simul_error = TRUE
  )
  
  target_missing_prop <- theta$prop_missing_data # 0.2
  tolerance <- 0.05 # 0.05 = between 15-25%
  
  for (g in seq_along(true_dataset$obs_dat)) {
    dat <- true_dataset$obs_dat[[g]]
    for (j in seq_len(ncol(dat))) {
      na_prop <- mean(is.na(dat[, j]))
      
      # Expect greater than 0.2 - tolerance
      expect_gt(
        na_prop,
        target_missing_prop - tolerance,
        label = paste("Group", g, "column", j, "has too few NAs:", na_prop)
      )
      
      # Expect less than 0.2 + tolerance
      expect_lt(
        na_prop,
        target_missing_prop + tolerance,
        label = paste("Group", g, "column", j, "has too many NAs:", na_prop)
      )
    }
  }
})

# ----------------------------------------------------------------------------

test_that("correct zeta is estimated", {
  set.seed(1)
  
  n_per_group <- rep(1000, n_groups)
  
  true_dataset <- simul_true_data(
    theta,
    n_per_group,
    range_dates,
    index_dates,
    simul_error = TRUE
  )
  
  expected_zeta <- theta$zeta
  tolerance <- 0.02
  
  for (g in seq_along(true_dataset$E)) {
    E_group <- true_dataset$E[[g]]
    
    for (j in seq_len(ncol(E_group))) {
      non_missing <- sum(E_group[, j] != -1)
      n_error <- sum(E_group[, j] == 1)
      
      zeta_est <- n_error / non_missing
        
        expect_true(
          abs(zeta_est - expected_zeta) < tolerance,
          paste("Group", g, "Column", j, 
                ": estimated zeta =", round(zeta_est, 3),
                "expected =", expected_zeta)
        )
      }
    }
})





