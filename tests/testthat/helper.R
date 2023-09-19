simulate_single_group <- function() {
  n_groups <- 1
  n_dates <- 2
  theta <- list()
  theta$prop_missing_data <- 0.05 ### this is missing from the estimation model (directly available from the data so not explicitely modelled)
  theta$zeta <- 0.05 ### probability that, when not missing, the date is recorded with error
  theta$mu <- list(5)
  theta$CV <- list(0.5)
  n_per_group <- 100
  range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
  index_dates <- list(matrix(c(1, 2), nrow = 2))

  siml <- simul_true_data(
    theta, n_per_group, range_dates, index_dates,
    simul_error = TRUE,
    remove_allNA_indiv = TRUE
  )
  all_NAs <- apply(siml$obs_dat[[1]], 1, function(row) all(is.na(row)))
  siml$obs_dat[[1]] <- siml$obs_dat[[1]][!all_NAs, ]
  list(obs_dat = siml$obs_dat, index_dates = index_dates)
}

simulate_four_groups <- function() {

    n_groups <- 4
    n_dates <- c(2, 3, 4, 4)
    theta <- list()
    theta$prop_missing_data <- 0.05 ### this is missing from the estimation model (directly available from the data so not explicitely modelled)
    theta$zeta <- 0.05 ### probability that, when not missing, the date is recorded with error
    theta$mu <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
    theta$CV <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))
    n_per_group <- rep(100, n_groups)
    range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
    index_dates <- list(matrix(c(1, 2), nrow = 2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
    index_dates_order <- list(
        matrix(c(1, 2), nrow = 2),
        cbind(c(1, 2), c(1, 3)),
        cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)),
        cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4))
    )
    siml <- simul_true_data(
        theta, n_per_group, range_dates, index_dates,
        simul_error = TRUE,
        remove_allNA_indiv = TRUE
    )
    all_NAs <- apply(siml$obs_dat[[1]], 1, function(row) all(is.na(row)))
    siml$obs_dat[[1]] <- siml$obs_dat[[1]][!all_NAs, ]
    list(obs_dat = siml$obs_dat, index_dates = index_dates)

}


bad_settings_one_group <- function() {
  mcmc_settings <- list(
    moves_switch = list(
      D_on = TRUE, E_on = TRUE,
      swapE_on = TRUE, mu_on = TRUE,
      CV_on = TRUE, zeta_on = TRUE
    ),
    moves_options = list(
        fraction_Di_to_update = 1/10,
        move_D_by_groups_of_size = 1,
        fraction_Ei_to_update = 1/10,
        sdlog_mu = list(0.05),
      sdlog_CV = list(0.25)
    ),
    init_options = list(mindelay = 0, maxdelay = 100),
    mcmc_settings = list(
        chain_properties = list(n_iter = 50, burnin = 10, record_every = 10)
    )
  )


  hyperparams <- list(
    shape1_prob_error = 3,
    shape2_prob_error = 12,
    mean_mean_delay = 100,
    mean_CV_delay = 100
  )

  list(mcmc_settings = mcmc_settings, hyperparameters = hyperparams)

}

## ngroups is either 1 or 4.
good_settings  <- function(ngroups = 1L) {
    if (ngroups == 1L) {
        settings <- bad_settings_one_group()
    } else {
        settings <- bad_settings_four_groups()
    }

  settings$mcmc_settings$chain_properties$n_iter <- 10000
  settings$mcmc_settings$chain_properties$burnin <- 5000
  settings
}


bad_settings_four_groups <- function() {
    list(
        moves_switch=list(D_on = TRUE, E_on = TRUE,  swapE_on = TRUE,  mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
        moves_options=list(
            fraction_Di_to_update = 1/10, move_D_by_groups_of_size = 1, fraction_Ei_to_update = 1/10,
            sdlog_mu = list(0.05, c(0.15, 0.15), c(0.15, 0.15, 0.15), c(0.25, 0.25, 0.25)),
            sdlog_CV = list(0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))),
         init_options = list(mindelay = 0, maxdelay = 100),
        mcmc_settings = list(
            chain_properties = list(
                n_iter = 500, burnin = 100, record_every = 10
            )
        )
    )
}


