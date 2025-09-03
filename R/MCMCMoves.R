###############################################
# Move functions
###############################################

#' Metropolis-Hastings accept/reject step
#' 
#' @description
#' Helper function to perform accept/reject step of a Metropolis-Hastings move.
#' It compares the final log posterior acceptance probability to a random draw
#' from a uniform distribution to decide whether to accept the proposed state
#' or keep the current state.
#' 
#' @param log_p_accept Log acceptance probability. Sum of the log posterior
#'  ratio and the log proposal correction factor. The function handles `NA`
#'  values by treating them as `-Inf` and caps positive values at `0`.
#' @param proposed_dat Augmented data representing the new, proposed state.
#'  This will be returned if the move is accepted.
#' @param current_dat Original state before the proposal. This will be returned
#'  if the move is rejected.
#'  
#' @return A list containing two elements:
#' \itemise{
#'  \item{`new_aug_dat`}: The augmented data set for the next step in the
#'     MCMC chain (either `proposed_dat` or `current_dat`).
#'   \item{`accept`}: An indicator with value `1` if the move was accepted or
#'     `0` if it was rejected.
#' }
#' 
#' @export
#' 
#' @examples
#' current_state <- list(D = 10)
#' proposed_state <- list(D = 15)
#' 
#' set.seed(1)
#' # Expect to accept
#' decide_acceptance(log_p_accept = 0, proposed_state, current_state)
#' 
#' # Expect to reject
#' decide_acceptance(log_p_accept = -10, proposed_state, current_state)
#' 
decide_acceptance <- function(log_p_accept, proposed_dat, current_dat) {
  if (is.na(log_p_accept)) log_p_accept <- -Inf
  if (log_p_accept > 0) log_p_accept <- 0
  
  if (log(runif(1)) < log_p_accept) {
    return(list(new_aug_dat = proposed_dat, accept = 1))
  } else {
    return(list(new_aug_dat = current_dat, accept = 0))
  }
}

# ----------------------------------------------------------------------------
# Move augmented dates D
# ----------------------------------------------------------------------------

#' Performs one iteration of an MCMC move for the augmented data
#'
#' @param i Index of individual(s) for whom augmented data should be moved.
#' @param group_idx Index of the group for whom augmented data should be moved.
#' @param date_idx Index of the date which should be moved.
#' @param curr_aug_dat The current augmented data; a list of observed data, in
#'  the format returned by \code{\link{simul_true_data}}.
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#'
#' @details \code{theta} should be a list containing:
#' \itemise{
#'  \item{\code{mu}: A list of length \code{n_groups} (the number of groups
#'   to be simulated data). Each element of \code{mu} should be a scalar of
#'   vector giving the mean delay(s) to use for simulation of dates in that
#'   group.}
#'  \item{\code{CV}: A list of length \code{n_groups}. Each element of
#'  \code{CV} should be a scalar of vector giving the coefficient o variation
#'   of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
#' \itemise{
#'  \item{\code{shape1_prob_error}: A scalar giving the first shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}: A scalar giving the second shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$CV}}
#' }
#' \code{index_dates} should be a list of length
#' \code{n_groups = length(obs_dat)}. Each element of \code{index_dates} should
#'  be a matrix with 2 rows and a number of columns corresponding to the delays
#'  of interest for that group. For each column (i.e. each delay), the first row
#'  gives the index of the origin date, and the second row gives the index of
#'  the destination date.
#'  The number of columns of index_dates[[k]] should match the length of
#'  theta$mu[[k]] and theta$CV[[k]].
#'
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and
#'  c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are
#'   respectively the mean and coefficient of variation of two delays: the
#'    first delay being between date 1 and date 2, and the second being between
#'     date 1 and date 3.
#'
#' The function performs the move as follows, using a Metropolis algorithm.
#' For the date to be moved, a new value is drawn from the marginal posterior
#' of one of the delays this date is involved in.
#' If the date is involved in several delays, one of the delays is randomly
#' selected.
#' The element E indicating whether the observed date is missing, recorded
#' correctly or recorded with error, is adjusted accordingly given the proposed
#' value of D.
#' The new augmented data is then accepted with probability given by the ratio
#' of the posterior values at the new augmented data and the old augmented data.
#'
#' @return A list of two elements:
#'  \itemise{
#'  \item{\code{new_aug_dat}: Same as \code{curr_aug_dat} but where the
#'  relevant dates have been updated}
#'  \item{\code{accept}: A scalar with value 1 if the move was accepted and
#'  0 otherwise}
#' }
#'
#' @export
#'
#' @examples
#' #'
#' # DONT THINK THESE ARE WORKING CORRECTLY
#' #Parameters
#' n_groups <- 4
#' n_per_group <- rep(100, n_groups)
#' n_dates <- c(2, 3, 4, 4)
#'
#' mu <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
#' cv <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))
#'
#' theta <- list(
#'   prop_missing_data = 0.2,
#'   zeta = 0.05,
#'   mu = mu,
#'   CV = cv
#' )
#'
#' hyperparameters <- list(
#'   # scalars giving the 1st and 2nd shape parameters for the beta prior for zeta
#'   shape1_prob_error = 3,
#'   shape2_prob_error = 12,
#'   # scalars giving the mean of the exponential prior used for mu and CV
#'   mean_mean_delay = 100,
#'   mean_CV_delay = 100
#' )
#'
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
#'                              as.Date("01/01/2015", "%d/%m/%Y")))
#'
#' index_dates <- list(
#'   matrix(c(1, 2), nrow = 2),
#'   cbind(c(1, 2), c(1, 3)),
#'   cbind(c(1, 2), c(2, 3),c(1, 4)),
#'   cbind(c(1, 2), c(2, 3), c(1, 4))
#' )
#'
#' # Simulate data
#' set.seed(1)
#' sim_data <- simul_true_data(theta,
#'                             n_per_group,
#'                             range_dates,
#'                             index_dates,
#'                             simul_error = TRUE,
#'                             remove_allNA_indiv = TRUE)
#'
#' obs_dat <- sim_data$obs_dat
#' mcmc_settings <- list(
#'   # moves_switch: booleans stating whether each parameter/augmented data
#'   # should be moved in the procedure or not.
#'   moves_switch = list(
#'     D_on = TRUE, # augmented dates (latent true dates)
#'     E_on = TRUE, # error indicators (-1, 0, 1)
#'     swapE_on = TRUE, # swaps error indicators - explore alternative errors
#'     mu_on = TRUE, # mean of each delay distribution
#'     CV_on = TRUE, # cv of each delay
#'     zeta_on = TRUE # probability of error
#'   ),
#'   # moves_options:
#'   moves_options = list(
#'     # Fraction of augmented dates to be updated at each iteration of the MCMC.
#'     fraction_Di_to_update = 1 / 10,
#'     # Number of augmented dates to be updated simultaneously in each group.
#'     move_D_by_groups_of_size = 1,
#'     # Fraction of indicators of whether observed dates are erroneous to be
#'     # updated at each iteration of the MCMC.
#'     fraction_Ei_to_update = 1 / 10,
#'     # List of SDs used for proposing moves of the mean delays of length n_groups.
#'     # Each element in the list should be a vector with length given by the
#'     # numbers of delays to be considered in this group.
#'     sdlog_mu = list(
#'       0.05,
#'       c(0.15, 0.15),
#'       c(0.15, 0.15, 0.15),
#'       c(0.25, 0.25, 0.25)
#'     ),
#'     # Same as above but for proposing moves of the CV of delays.
#'     sdlog_CV = list(
#'       0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
#'   ),
#'   # minimum and maximum delays, below/above which dates are considered
#'   # incompatible with one another at the initialisation stage of the MCMC.
#'   init_options = list(mindelay = 0, maxdelay = 100),
#'   # total number of iterations, initial burnin and then after burnin how many
#'   # iterations should be recorded (thinning). (500 - 50) / 10 = 45 samples from
#'   # the posterior for each dataset.
#'   chain_properties = list(n_iter = 500, burnin = 50, record_every = 10)
#' )
#' curr_aug_dat <- initialise_aug_data(obs_dat, index_dates, MCMC_settings = mcmc_settings)
#' theta <- initialise_theta_from_aug_dat(curr_aug_dat, index_dates)
#'
#' # Example where date missing (E = -1)
#' group_idx <- 1
#' i <- 8
#' date_idx <- 2
#'
#' curr_aug_dat$D[[group_idx]][i, date_idx]
#' curr_aug_dat$E[[group_idx]][i, date_idx]
#'
#' # Move a date for individual 8, group 1, date index 1
#' set.seed(1)
#' result1 <- move_Di(i, group_idx, date_idx,
#'                   curr_aug_dat = curr_aug_dat,
#'                   theta = theta,
#'                   obs_dat = obs_dat,
#'                   hyperparameters = hyperparameters,
#'                   index_dates = index_dates)
#'
#' # Check result
#' result1$accept                                  # 1 = accepted, 0 = rejected
#' result1$new_aug_dat$D[[group_idx]][i, date_idx] # New proposed date value
#' curr_aug_dat$D[[group_idx]][i, date_idx]       # Old proposed date value
#' obs_dat[[group_idx]][i, date_idx]              # Original observed date missing
#'
#' # Example where date observed with error (E = 1)
#' obs_dat <- sim_data$obs_dat
#' curr_aug_dat <- initialise_aug_data(obs_dat, index_dates, MCMC_settings = mcmc_settings)
#' theta <- initialise_theta_from_aug_dat(curr_aug_dat, index_dates)
#' group_idx <- 1
#' i <- 20
#' date_idx <- 1
#'
#' curr_aug_dat$E[[group_idx]][i, date_idx]
#'
#' # Move a date for individual 20, group 1, date index 1
#' set.seed(10)
#' result2 <- move_Di(i, group_idx, date_idx,
#'                   curr_aug_dat = curr_aug_dat,
#'                   theta = theta,
#'                   obs_dat = obs_dat,
#'                   hyperparameters = hyperparameters,
#'                   index_dates = index_dates)
#'
#' # Check result
#' result2$accept                                  # 1 = accepted, 0 = rejected
#' result2$new_aug_dat$D[[group_idx]][i, date_idx] # New proposed date value
#' curr_aug_dat$D[[group_idx]][i, date_idx]       # Old proposed date value
#' obs_dat[[group_idx]][i, date_idx]              # Original observed date
#' #'
move_Di <- function(i,
                    group_idx,
                    date_idx,
                    curr_aug_dat,
                    theta,
                    obs_dat,
                    hyperparameters,
                    index_dates,
                    range_dates = NULL) {

  if (is.null(range_dates)) range_dates <- find_range(obs_dat)

  tmp <- propose_new_delay(i,
                           group_idx,
                           date_idx,
                           curr_aug_dat,
                           theta,
                           obs_dat,
                           hyperparameters,
                           index_dates,
                           range_dates)

  proposed_aug_dat <- tmp$proposed_aug_dat
  curr_delay <- tmp$curr_delay
  sample_delay <- tmp$sample_delay

  # Calculate posterior ratio -------------------------------------------------
  # probability of acceptance = log P(new) - log P(old)

  # index for delays that are affected by the change in date date_idx
  delay_idx <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)[, 2]

  ll_proposed <- LL_observation_term_by_group_delay_and_indiv(
    proposed_aug_dat, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  )

  ll_current <- LL_observation_term_by_group_delay_and_indiv(
    curr_aug_dat, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  )

  ratio_post <- ll_proposed - ll_current

  # Add error term difference only if E changed
  different_E <- proposed_aug_dat$E[[group_idx]][i, date_idx] !=
    curr_aug_dat$E[[group_idx]][i, date_idx]

  if (any(different_E)) {
    ratio_post <- ratio_post +
      LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta,
                                             group_idx, date_idx, i) -
      LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta,
                                             group_idx, date_idx, i)
  }

  # Add delay likelihood differences for each affected delay
  for (d in delay_idx) {
    ratio_post <- ratio_post +
      LL_delays_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat,
                                              group_idx, d, i, index_dates) -
      LL_delays_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat,
                                              group_idx, d, i, index_dates)
  }

  ratio_post <- sum(ratio_post)

  ### note that ratio_post should be the same as:
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat,
  # hyperparameters, index_dates, range_dates) -
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates, range_dates)

  # Proposal correction factor ------------------------------------------------

  # Correction factor needed as this move is not symmetrical
  # where Q = proposal distribution, theta_old = curr_delay and theta_new = sample_delay:
  # corr = Q(theta_old | theta_new) / Q(theta_new | theta_old)
  # log_corr = log(Q(theta_old | theta_new)) - log(Q(theta_new | theta_old))

  tmp <- get_correct_factor_new_delay(curr_delay, sample_delay, theta, group_idx, tmp$which_delay)
  ratio_prop <- tmp["prob_proposing_curr_value"] - tmp["prob_proposing_new_value"]

  # Acceptance probability ----------------------------------------------------
  p_accept <- ratio_post + ratio_prop

  # Accept/reject the proposal ------------------------------------------------
  res <- decide_acceptance(log_p_accept = p_accept,
                           proposed_dat = proposed_aug_dat,
                           current_dat = curr_aug_dat)
  
  return(res)

}
# test_move_Di <- move_Di(i=1, group_idx=1, date_idx=1, curr_aug_dat = aug_dat,
# theta, obs_dat, hyperparameters)
# test_move_Di$new_aug_dat$D[[1]][1,1] # new value
# aug_dat$D[[1]][1,1] # old value


# -----------------------------------------------------------------------------
# Move augmented indicator for whether date is correctly recorded E
# -----------------------------------------------------------------------------

#' Propose a new date when moving from no error (`E=0`) to observed with error
#'  (`E=1`), ensuring the new date is different from the observed date.
#'  
#' @description Helper function that proposes a new true date `D` when an error
#' indicator `E` is being swapped from 0 to 1. The proposed date is sampled
#' from one of the delay distributions the date is involved in.
#'
#' @details The function first identifies all delays that the specified date
#' (`date_idx`) is a part of. It then chooses one of these delays at random.
#' A new delay value is sampled from the corresponding discrete gamma
#' distribution, and this is used to calculate a new proposed date `D`. The
#' function ensures that the proposed date is not the same as the observed date,
#' which is necessary for the date to have an error indicator of `E=1`. This
#' function only proposes a new `D` and doesn't update `E`.
#'
#' @param i Index of the individual for whom augmented data should be moved.
#' @param group_idx Index of the group the individual is in.
#' @param date_idx Numeric vector of date indices being moved.
#' @param curr_aug_dat The current augmented data.
#' @param theta List of parameters, including `mu`, `CV` and `zeta`.
#' @param obs_dat List of observed data.
#' @param hyperparameters List of hyperparameters.
#' @param index_dates List containing the delays defined for each group.
#' @param range_dates Vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#'
#' @return Numeric vector of proposed new dates `D`.
#' @export
#' 
#' @examples
#' # Create minimal data for a single individual
#' theta_ex <- list(mu = list(c(5, 10)), CV = list(c(0.2, 0.3)))
#' index_dates_ex <- list(matrix(c(1, 2, 2, 3), nrow = 2, byrow = TRUE))
#' curr_aug_dat_ex <- list(D = list(matrix(c(10, 15, 25), nrow = 1)))
#' obs_dat_ex <- list(matrix(c(10, 16, 25), nrow = 1))
#' 
#' # Set seed for reproducible random choice of delay
#' set.seed(1)
#' 
#' # Propose a new value for the second date (date_idx = 2)
#' # This date is involved in two delays. With seed(1), the first delay (1->2)
#' # is chosen. The new date will be date 1 + sample_delay = 10 + 4 = 14.
#' new_date <- propose_move_from_E0_to_E1(
#'   i = 1, group_idx = 1, date_idx = 2,
#'   curr_aug_dat = curr_aug_dat_ex,
#'   theta = theta_ex,
#'   obs_dat = obs_dat_ex,
#'   hyperparameters = NULL,
#'   index_dates = index_dates_ex
#' )
#' print(new_date)
#' 
propose_move_from_E0_to_E1 <- function(i,
                                       group_idx,
                                       date_idx,
                                       curr_aug_dat,
                                       theta,
                                       obs_dat,
                                       hyperparameters,
                                       index_dates,
                                       range_dates = NULL) {

  # Find all delays these dates are involved in
  x <- lapply(seq_along(date_idx), function(e) {
    which(index_dates[[group_idx]] == date_idx[e], arr.ind = TRUE)
  })

  which_delay <- lapply(seq_along(date_idx), function(e) x[[e]][, 2])

  # ANNE: what are the paired dates involved in those delays
  # ANNE: what are their index
  from_idx <- lapply(seq_along(date_idx), function(e) {
    x_e <- x[[e]]
    sapply(seq_len(nrow(x_e)), function(k) {
      index_dates[[group_idx]][-x_e[k, 1], x_e[k, 2]]
    })
  })

  # ANNE: what are the actual dates
  # use from_idx to extract corresponding D
  from_value <- lapply(seq_along(date_idx), function(e) {
    x_e <- x[[e]]
    sapply(seq_len(nrow(x_e)), function(k) {
      curr_aug_dat$D[[group_idx]][i, index_dates[[group_idx]][-x_e[k, 1], x_e[k, 2]]]
    })
  })

  # if several delays involved, choose one at random
  tmp <- lapply(seq_along(date_idx), function(e) {
    sample(seq_along(from_idx[[e]]), 1)
  })

  which_delay <- lapply(seq_along(date_idx), function(e) {
    which_delay[[e]][tmp[[e]]]
  })

  from_idx <- sapply(seq_along(date_idx), function(e) {
    from_idx[[e]][tmp[[e]]]
  })

  from_value <- sapply(seq_along(date_idx), function(e) {
    from_value[[e]][tmp[[e]]]
  })

  # ANNE: this function just samples from the selected delay to obtain a new
  # proposed augmented date D
  get_one_proposed_aug_value <- function(e) {
    mu <- theta$mu[[group_idx]][which_delay[[e]]]
    cv <- theta$CV[[group_idx]][which_delay[[e]]]

    sample_delay <- discr_gamma_sample(1, mu = mu, cv = cv)

    if (date_idx[e] < from_idx[e]) {
      proposed_aug_dat_value <- from_value[e] - sample_delay
    } else {
      proposed_aug_dat_value <- from_value[e] + sample_delay
    }

    # while we haven't moved to a situation where E = 1, try again
    while (proposed_aug_dat_value == obs_dat[[group_idx]][i, date_idx[e]]) {
      sample_delay <- discr_gamma_sample(1, mu = mu, cv = cv)

      if (date_idx[e] < from_idx[e]) {
        proposed_aug_dat_value <- from_value[e] - sample_delay
      } else {
        proposed_aug_dat_value <- from_value[e] + sample_delay
      }
    }

    proposed_aug_dat_value
  }

  proposed_aug_dat_value <- sapply(
    seq_len(length(date_idx)), get_one_proposed_aug_value
  )

  return(proposed_aug_dat_value)
}


#' Propose the observed date when moving from error (`E=1`) to no error (`E=0`).
#' 
#' @description Helper function that proposes a move to the originally observed
#'  date from `obs_dat` when `E` is being swapped from 1 to 0.
#'
#' @details This is a deterministic move. When an error is corrected, the new
#' "true" date is assumed to be the date that was originally observed.
#' This function operates on a single `date_idx` at a time. The calling
#' function (`perform_E1_to_E0_swap`) is responsible for looping
#' through multiple dates if necessary.
#'
#' @param i Index of the individual.
#' @param group_idx Index of the group.
#' @param date_idx Index of the date being moved.
#' @param obs_dat A list of observed data.
#'
#' @return Proposed true date, which is the same as the observed date.
propose_move_from_E1_to_E0 <- function(i,
                                       group_idx,
                                       date_idx,
                                       obs_dat) {

  proposed_aug_dat_value <- obs_dat[[group_idx]][i, date_idx]

  proposed_aug_dat_value
}


# ----------------------------------------------------------------------------
# Combine compute_p_accept_move_* functions
# ----------------------------------------------------------------------------

#' Computes acceptance probability for a move of the error indicator E
#'
#' @description Calculates the log acceptance probability for an MCMC
#' move that flips an error indicator from E=0 to E=1 or vice-versa. It
#' computes both the log posterior ratio and the proposal correction factor
#' for the asymmetric move.
#'
#' @param direction A character string, either "E0_to_E1" or "E1_to_E0",
#'   specifying the direction of the move.
#' @param i Index of the individual for whom the move is being made.
#' @param group_idx Index of the group for the move.
#' @param date_idx Index of the date for the move.
#' @param curr_aug_dat The current augmented data before the proposal.
#' @param proposed_aug_dat The proposed augmented data after the move.
#' @param theta List of model parameters.
#' @param obs_dat A list of observed data.
#' @param index_dates A list defining the delays for each group.
#' @param range_dates A vector containing the range of dates in `obs_dat`.
#'   If NULL, it will be computed automatically.
#'
#' @return A vector of length 2. The first element is the log posterior
#'   difference (proposed - current), and the second is the log proposal
#'   correction factor.
#'
compute_p_accept_move_E <- function(direction = c("E0_to_E1", "E1_to_E0"),
                                    i,
                                    group_idx,
                                    date_idx,
                                    curr_aug_dat,
                                    proposed_aug_dat,
                                    theta,
                                    obs_dat,
                                    index_dates,
                                    range_dates = NULL) {
  
  if (!direction %in% c("E0_to_E1", "E1_to_E0")) {
    stop("direction must be either `E0_to_E1` or `E1_to_E0`")
  }
  
  direction <- match.arg(direction)
  
  # Delays that are affected by changing this specific date (date_idx)
  delay_idx <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)[, 2]
  
  # Compute log posterior ratio ----------------------------------------------
  # Calculate the log(P(new) / P(old)) by summing the changes in the
  # log-likelihood for each component
  
  # Difference in observation likelihood
  ratio_post <- LL_observation_term_by_group_delay_and_indiv(
    proposed_aug_dat, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  ) - LL_observation_term_by_group_delay_and_indiv(
    curr_aug_dat, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  )
  
  # Difference in error likelihood
  ratio_post <- ratio_post + LL_error_term_by_group_delay_and_indiv(
    proposed_aug_dat, theta, group_idx, date_idx, i
  ) - LL_error_term_by_group_delay_and_indiv(
    curr_aug_dat, theta, group_idx, date_idx, i
  )
  
  # For each affected delay, add the difference in delay likelihood
  for (d in delay_idx) {
    ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(
      proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates
    ) - LL_delays_term_by_group_delay_and_indiv(
      curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates
    )
  }
  
  # Combine all log-likelihood differences
  ratio_post <- sum(ratio_post)
  
  ### note that ratio_post should be the same as:
  ## ANNE: TODO check this works
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat,
  # hyperparameters, index_dates) -
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates)
  
  # Proposal correction factor -----------------------------------------------
  # Calculate the log of the posterior ratio: log(q(old|new) / q(new|old))
  # to correct for the asymmetric proposal
  
  # Index for the date within index_dates
  x <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)
  which_delay <- x[, 2]
  
  # Index for the other date involved in each delay
  from_idx <- sapply(
    seq_len(nrow(x)), function(k) index_dates[[group_idx]][-x[k, 1], x[k, 2]]
  )
  
  # Extract value of the other date for each delay
  from_value <- sapply(
    seq_len(nrow(x)),
    function(k) curr_aug_dat$D[[group_idx]][i, index_dates[[group_idx]][-x[k, 1], x[k, 2]]]
  )
  
  # To calculate the correction, take the date value from the state where E=1
  date_value_for_corr <- if (direction == "E0_to_E1") {
    # In this move, the proposed state has E=1 so we use this date
    proposed_aug_dat$D[[group_idx]][i, date_idx]
  } else if (direction == "E1_to_E0") {
    # In this move, the current state has E=1 so we use this date
    curr_aug_dat$D[[group_idx]][i, date_idx]
  }
  
  # Calculate the probability (K) of proposing the 'date_value_for_corr' for
  # a single delay it is involved in
  find_correction_factor <- function(e) {
    
    # Calculate delay ensuring it's a positive value
    if (date_idx < from_idx[e]) {
      delay <- from_value[e] - date_value_for_corr
    } else {
      delay <- date_value_for_corr - from_value[e]
    }
    
    # Delay corresponding to the observed date. Moving to E=1 requires
    # D != obs_dat so this specific delay value is 'forbidden'
    if (date_idx < from_idx[e]) {
      forbidden_delay <- from_value[e] - obs_dat[[group_idx]][i, date_idx]
    } else {
      forbidden_delay <- obs_dat[[group_idx]][i, date_idx] - from_value[e]
    }
    
    mu <- theta$mu[[group_idx]][which_delay[e]]
    cv <- theta$CV[[group_idx]][which_delay[e]]
    
    # Calculate the probability of randomly drawing the proposed delay,
    # renormalising the distribution because we cannot choose forbidden delay
    K <- DiscrGamma(k = delay, mu = mu, cv = cv, log = FALSE) /
      (1 - DiscrGamma(k = forbidden_delay, mu = mu, cv = cv, log = FALSE))
    
    return(K)
  }
  
  # Average correction factors over all delays involving this date
  # ANNE: accounts for the fact that the proposal distrib is a mixture distrib
  # where the index of the delay used for sampling is drawn at random.
  # hence the correction factor, which should be the probability of drawing the
  # proposed delay irrespective of which delay index was chosen
  # hence it's P(proposed delay | delay index 1) * P(delay index 1) +
  # P(proposed delay | delay index 2) * P(delay index 2) + ...
  # which simplifies to just the mean of P(proposed delay | delay index i)
  # because all the P(delay index i) are the same
  K <- mean(sapply(seq_along(which_delay), find_correction_factor))
  
  # Log proposal correction factor - calculated differently depending on the
  # direction of the move, because the proposal mechanism is asymmetric.
  #
  # E1 -> E0: There is only 1 way to move (true date must become observed date).
  # The probability of proposing this specific move is 1, so it's
  # log probability is 0.
  #
  # E0 -> E1: The probability of proposing a specific new date is calculated
  # from the delay distribution (DiscrGamma), discounting the one value that
  # would correspond to E=0.
  #
  # The log correction is always logP(proposing current) - logP(proposing new)
  logcorrection <- if (direction == "E0_to_E1") {
    # Forward move is probabilistic (K), reverse is deterministic (1)
    # log(1) - log(K) = -log(K)
    -log(K)
  } else if (direction == "E1_to_E0") {
    # Forward move is deterministic (1), reverse is probabilistic (K)
    # log(K) - log(1) = +log(K)
    +log(K)
  }

  return(c(ratio_post, logcorrection))
}

# ----------------------------------------------------------------------------
# Move E
# ----------------------------------------------------------------------------

#' Performs one iteration of an MCMC move for the augmented data representing
#'  the indicator of error in observations
#'
#' @param i The index of the individual for whom augmented data should be moved
#' @param group_idx The index of the group for whom augmented data should be
#'  moved
#' @param date_idx The index of the date for which the error indicator which
#'  should be moved
#' @param curr_aug_dat The current augmented data; a list of observed data, in
#'  the format returned by \code{\link{simul_true_data}}.
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#' @details \code{theta} should be a list containing:
#' \itemise{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups
#'   to be simulated data). Each element of \code{mu} should be a scalar of
#'    vector giving the mean delay(s) to use for simulation of dates in that
#'     group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar of vector giving the coefficient o variation
#'    of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
#' \itemise{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape
#'   parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$CV}}
#' }
#' \code{index_dates} should be a list of length
#'  \code{n_groups = length(obs_dat)}. Each element of \code{index_dates}
#'   should be a matrix with 2 rows and a number of columns corresponding to
#'    the delays of interest for that group. For each column (i.e. each delay),
#'     the first row gives the index of the origin date, and the second row
#'      gives the index of the destination date.
#' The number of columns of index_dates[[k]] should match the length of
#'  theta$mu[[k]] and theta$CV[[k]]
#'
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and
#'  c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are
#'   respectively the mean and coefficient of variation of two delays: the
#'    first delay being between date 1 and date 2, and the second being between
#'     date 1 and date 3.
#'
#' The function performs the move as follows, using a Metropolis-Hastings
#'  algorithm.
#' If E=-1 nothing happens.
#' If E=1, we propose a move to E=0 and hence D=the observed data.
#' If E=0, we propose a move to E=1. D is then moved as follows: a new value is
#'  drawn from the marginal posterior of one of the delays this date is
#'   involved in, repeatedly until D falls on a different day than the observed
#'    date to be consistent with E=1.
#' If the date is involved in several delays, one of the delays is randomly
#'  selected.
#' This move is not symmetrical so we use a correction factor in computing the
#'  probability of acceptance in the Metropolis Hastings which accounts for the
#'   asymetry.
#' @return A list of two elements:
#'  \itemise{
#'  \item{\code{new_aug_dat}}{: Same as \code{curr_aug_dat} but where the
#'   relevant indicators of errors in dates have been updated}
#'  \item{\code{accept}}{: A scalar with value 1 if the move was accepted and
#'   0 otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
move_Ei <- function(i,
                    group_idx,
                    date_idx,
                    curr_aug_dat,
                    theta,
                    obs_dat,
                    hyperparameters,
                    index_dates,
                    range_dates = NULL) {
  
  if (length(i) > 1) {
    i <- i[1]
    warning("In move_Ei, i should be a single number, not a vector. Using i[1] instead.")
  }

  if (is.null(range_dates)) range_dates <- find_range(obs_dat)

  curr_E_value <- curr_aug_dat$E[[group_idx]][i, date_idx]

  # if date is missing do nothing
  if (curr_E_value == -1) {
    return(list(new_aug_dat = curr_aug_dat, accept = 0))
  }
  
  proposed_aug_dat <- curr_aug_dat
  p_accept <- -Inf # reject as the default
  
  # moving from E=0 to E=1
  if (curr_E_value == 0) {
  
    proposed_aug_dat$E[[group_idx]][i, date_idx] <- 1
    proposed_aug_dat$D[[group_idx]][i, date_idx] <- propose_move_from_E0_to_E1(
      i, group_idx, date_idx, curr_aug_dat, theta,
      obs_dat, hyperparameters, index_dates, range_dates
    )

    tmp <- compute_p_accept_move_E(direction = "E0_to_E1",
      i, group_idx, date_idx, curr_aug_dat, proposed_aug_dat,
      theta, obs_dat, index_dates, range_dates
    )

    if (!any(is.infinite(tmp))) p_accept <- sum(tmp)
  
    # moving from E=1 to E=0
  } else if (curr_E_value == 1) {
    
    proposed_aug_dat$E[[group_idx]][i, date_idx] <- 0
    proposed_aug_dat$D[[group_idx]][i, date_idx] <- propose_move_from_E1_to_E0(
      i, group_idx, date_idx, obs_dat
    )
    
    tmp <- compute_p_accept_move_E(direction = "E1_to_E0",
      i, group_idx, date_idx, curr_aug_dat, proposed_aug_dat,
      theta, obs_dat, index_dates, range_dates
    )
    
    if (!any(is.infinite(tmp))) p_accept <- sum(tmp)
  }

    res <- decide_acceptance(p_accept, proposed_aug_dat, curr_aug_dat)
    
    # return a list of size 2 where
    #		the first value is the new augmented data set in the chain
    #		the second value is 1 if the proposed value was accepted, 0 otherwise
    return(res)
    
  }

# ----------------------------------------------------------------------------
# Swap Es - when only 2 Es, related by a delay, are recorded, one with error
# and one without error, propose to swap the two
# ----------------------------------------------------------------------------

#' Find individuals with mixed error types for swapping
#' 
#' @description
#' Helper function for the MCMC sampler that identifies which individuals are
#' eligible for the `swap_Ei` move. An individual is eligible if their set of
#' recorded dates contains at least one date marked as correct (`E=0`) and at
#' least one date marked as an error (`E=1`). Individuals with only correct
#' dates, only erroneous dates, or only missing dates, are not eligible for
#' this swap move.
#' 
#' @param group_idx Integer specifying the index of the group to check
#' @param curr_aug_dat Current augmented data list, which must contain a list
#'  of error matrices named `E`
#'  
#' @return Vector containing the row indices of the individuals (within the
#'  specified group) who have a mix of `E=0` and `E=1` values. Returns an
#'  empty vector if no individuals are found.
#' 
#' @export
#' 
#' @examples
#' # Augmented data with one group and four individuals
#' aug_dat <- list(E = list(
#'   matrix(c(
#'     0, 0, -1, 0,  # all correct or missing (not eligible)
#'     1, 1, 1, -1,   # all errors or missing (not eligible)
#'     0, 1, -1, 0,  # mix of 0, 1, and -1 (eligible)
#'     1, 0, 1, 0   # mix of 0 and 1 (eligible)
#'   ), nrow = 4, byrow = TRUE)
#' ))
#'
#' # find_Eis_to_swap should identify individuals 3 and 4 as eligible
#' find_Eis_to_swap(group_idx = 1, curr_aug_dat = aug_dat)
#' 
find_Eis_to_swap <- function(group_idx, curr_aug_dat) {
  
  Es <- curr_aug_dat$E[[group_idx]]
  
  # for each row (individual) check for more than one unique non_missing entry
  which(sapply(seq_len(nrow(Es)), function(i) {
    non_missing <- Es[i, Es[i, ] != -1]
    unique_values_of_E_i <- unique(non_missing)
    # length only >1 if both a 0 and 1 present
    length(unique_values_of_E_i) > 1
  }))
}

## FORWARD MOVES ------------------------------------------------------------

#' Perform the E=1 to E=0 swap for a set of dates
#' 
#' @description
#' First sequential step within `swap_Ei`. It takes all dates for a given
#' individual that are currently marked as an error (`E=1`) and proposes a new
#' state where they are marked as correct (`E=0`).
#' 
#' @details
#' This is a deterministic proposal. For each date being moved, the error
#' indicator `E` is set to `0`, and the true date `D` is set to the
#' corresponding observed date from `obs_dat`. As this move is deterministic,
#' the forward proposal probability is 1, and the log-probability is 0.
#' 
#' @param i Index of the individual.
#' @param group_idx Index of the group.
#' @param date_idx_E1_to_E0 Numeric vector of column indices for the dates
#'   that are `E=1` and will be moved.
#' @param curr_aug_dat Current augmented data list before the move.
#' @param obs_dat List of the observed data.
#' @param theta List of model parameters (not used in this specific function
#'   but passed for consistency with other swap steps).
#' @param hyperparameters List of model hyperparameters (not used in this
#'   specific function).
#' @param index_dates List defining delays for each group (not used in this
#'   specific function).
#' @param range_dates Vector of the overall date range (not used in this
#'   specific function).
#'
#' @return The augmented data list (`proposed_aug_dat`) with the relevant dates
#'   and error indicators updated.
#' @export
#' @seealso `swap_Ei`, `propose_move_from_E1_to_E0`
#'
perform_E1_to_E0_swap <- function(i, group_idx, date_idx_E1_to_E0,
                                  curr_aug_dat, theta, obs_dat,
                                  hyperparameters, index_dates, range_dates) {
  
  proposed_aug_dat <- curr_aug_dat
  if (length(date_idx_E1_to_E0) == 0) {
    return(proposed_aug_dat)
  }
  
  for (k in date_idx_E1_to_E0) {
    proposed_aug_dat$E[[group_idx]][i, k] <- 0
    proposed_aug_dat$D[[group_idx]][i, k] <-
      propose_move_from_E1_to_E0(i, group_idx, k, obs_dat)
  }
  return(proposed_aug_dat)
}

#' Perform the E=0 to E=1 swap for a set of dates
#' 
#' @description
#' Second sequential step within `swap_Ei`. It takes all dates for a given
#' individual that are currently marked as correct (`E=0`) and proposes a new
#' state where they are marked as an error (`E=1`).
#' 
#' @details
#' This is a non-deterministic proposal. For each date being moved, the error
#' indicator `E` is set to `1`, and a new plausible true date `D` is sampled
#' via `propose_move_from_E0_to_E1`. This sampling is conditional on the other
#' dates, including those updated in the previous step of the swap.
#' 
#' This function also calculates and accumulates the log-probability
#' of making these sequential proposals. This is done by calling
#' `compute_p_accept_move_E` for each move and summing the resulting
#' correction factors.
#' 
#' @param i Index of the individual.
#' @param group_idx Index of the group.
#' @param date_idx_E0_to_E1 Numeric vector of column indices for the dates
#'   that are `E=0` and will be moved.
#' @param current_aug_dat Augmented data list after the E=1 to E=0 swap
#'   has already been performed in the preceding step.
#' @param theta List of model parameters (`mu`, `CV`, `zeta`).
#' @param obs_dat List of the observed data.
#' @param hyperparameters List of the model's hyperparameters.
#' @param index_dates List defining the delays for each group.
#' @param range_dates Vector containing the overall date range.
#' 
#' @return A list containing two elements:
#' \itemise{
#'   \item{`proposed_aug_dat`}: Augmented data with the relevant dates
#'     updated to the new `E=1` state.
#'   \item{`correction_factor`}: Sum of the log-probabilities for the
#'     proposal of these dates.
#' }
#' @export
#' @seealso `swap_Ei`, `propose_move_from_E0_to_E1`, `compute_p_accept_move_E`
#'
perform_E0_to_E1_swap <- function(i,
                                  group_idx,
                                  date_idx_E0_to_E1,
                                  current_aug_dat,
                                  theta,
                                  obs_dat,
                                  hyperparameters,
                                  index_dates,
                                  range_dates) {
  
  proposed_aug_dat <- current_aug_dat
  correction_factor <- 0
  
  # If there are no dates to swap return inputs unchanged
  if (length(date_idx_E0_to_E1) == 0) {
    return(list(
      proposed_aug_dat = proposed_aug_dat,
      correction_factor = correction_factor
    ))
  }
  
  # Loop through each date that was originally E=0
  for (k in date_idx_E0_to_E1) {
    state_before_swap <- proposed_aug_dat
    
    # Propose new state for this date
    proposed_aug_dat$E[[group_idx]][i, k] <- 1
    proposed_aug_dat$D[[group_idx]][i, k] <-
      propose_move_from_E0_to_E1(
        i, group_idx, k, state_before_swap,
        theta, obs_dat, hyperparameters, index_dates, range_dates
      )
    
    # Calculate and accumulate correction factors
    correction_factor <- correction_factor +
      compute_p_accept_move_E(direction = "E0_to_E1",
        i = i,
        group_idx = group_idx,
        date_idx = k,
        curr_aug_dat = state_before_swap,
        proposed_aug_dat = proposed_aug_dat,
        theta = theta,
        obs_dat = obs_dat,
        index_dates = index_dates,
        range_dates = range_dates
      )[2] # only need second element of the result
  }
  
  return(list(
    proposed_aug_dat = proposed_aug_dat,
    correction_factor = correction_factor
  ))
}

#' Resample missing dates
#' 
#' @description
#' Resamples any missing dates (`E=-1`) within the `swap_Ei` function to ensure
#' they remain consistent with the newly proposed dates from the preceding
#' steps of the swap.
#' 
#' @details
#' After the error indicators (`E=0` and `E=1`) have been swapped, this function
#' iterates through dates marked as missing (`E=-1`) for a given individual. For
#' each missing date it calls `propose_new_delay()` to generate a new date that
#' is conditional on the already updated dates. It also calculates the log-
#' probability of making these proposals.
#' 
#' @param i Index of the individual for whom the dates are being resampled
#' @param group_idx Index of the group for the individual
#' @param date_idx_resample Numeric vector of column indixes for the missing
#'  dates (`E=-1`) that need to be resampled.
#' @param current_aug_dat Augmented data list after the E=0/E=1 swaps that have
#'  already been performed in the preceding steps.
#' @param theta List of model parameters (`mu`, `CV`, `zeta`)
#' @param obs_dat List of observed data
#' @param hyperparameters List of model hyperparameters
#' @param index_dates List defining the delays for each group
#' @param range_dates Vector containing the overall range of dates to consider
#' 
#' @return List containing two elements:
#' \itemise{
#'  \item{`proposed_aug_dat`}: Augmented data with missing dates resampled to
#'    be consistent with the other proposed changes.
#'  \item{`correction_factor`}: Sum of the log-probabilities for the proposal
#'    of the missing dates.
#' }
#' @export
#' 
resample_missing_dates <- function(i,
                                   group_idx,
                                   date_idx_resample,
                                   current_aug_dat,
                                   theta,
                                   obs_dat,
                                   hyperparameters,
                                   index_dates,
                                   range_dates) {
  
  proposed_aug_dat <- current_aug_dat
  correction_factor <- 0
  
  # If there are no dates to resample return inputs unchanged
  if (length(date_idx_resample) == 0) {
    return(list(
      proposed_aug_dat = proposed_aug_dat,
      correction_factor = correction_factor
    ))
  }
  
  # Loop through each date that was originally missing
  for (k in date_idx_resample) {
    # Propose a new delay and get the intermediate values
    tmp_delay_info <- propose_new_delay(
      i, group_idx, k, proposed_aug_dat,
      theta, obs_dat, hyperparameters, index_dates, range_dates
    )
    
    # Update the augmented data with the new proposal from this iteration
    proposed_aug_dat <- tmp_delay_info$proposed_aug_dat
    
    # Calculate and accumulate the correction factor for this specific move
    correction_factor <- correction_factor +
      get_correct_factor_new_delay(
        tmp_delay_info$curr_delay,
        tmp_delay_info$sample_delay,
        theta,
        group_idx,
        tmp_delay_info$which_delay
      )["prob_proposing_new_value"]
  }
  
  return(list(
    proposed_aug_dat = proposed_aug_dat,
    correction_factor = correction_factor
  ))
}


## REVERSE MOVES ------------------------------------------------------------

#' Calculate reverse proposal probability for resampled missing dates
#' 
#' @description
#' Sequential step within `swap_Ei`. It performs the reverse of the
#' `resample_missing_dates` step. It reverts the dates for missing indices
#' back to their original values and calculates the correction factor for this
#' reverse move.
#' 
#' @param i Index of the individual.
#' @param group_idx Index of the group.
#' @param date_idx_resample Numeric vector of column indices for the dates
#'   that were resampled in the forward move.
#' @param final_proposed_dat Augmented data list after the full forward
#'   proposal has been completed.
#' @param original_dat Augmented data list before any proposals were made.
#' @param theta List of model parameters (`mu`, `CV`, `zeta`).
#' @param index_dates List defining the delays for each group.
#' 
#' @return A list containing two elements:
#' \itemise{
#'   \item{`reverted_aug_dat`}: Augmented data with the resampled dates
#'     reverted back to their original state.
#'   \item{`correction_factor`}: Sum of the log-probabilities for the
#'     reverse proposal of the missing dates.
#' }
#' @export
#' @seealso `resample_missing_dates`
#' 
reverse_resample_missing_dates <- function(i,
                                           group_idx,
                                           date_idx_resample,
                                           final_proposed_dat,
                                           original_dat,
                                           theta,
                                           index_dates) {
  
  correction_factor_rev <- 0
  reverted_dat <- final_proposed_dat
  
  if (length(date_idx_resample) > 0) {
    # Delays from final proposed state
    curr_delays <- compute_delta(final_proposed_dat$D, index_dates)[[group_idx]][i, ]
    
    # Revert the dates for all resampled incides at once
    reverted_dat$D[[group_idx]][i, date_idx_resample] <-
      original_dat$D[[group_idx]][i, date_idx_resample]
    
    # Get new delays after reverting the dates
    new_delays <- compute_delta(reverted_dat$D, index_dates)[[group_idx]][i, ]
    
    # Calculate the reverse correction factor
    for (k in date_idx_resample) {
      # Find which delays were affected by changing date k
      tmp_delay_idx <- which(index_dates[[group_idx]] == k, arr.ind = TRUE)[, 2]
      for (kk in tmp_delay_idx) {
        correction_factor_rev <- correction_factor_rev +
          (1 / length(tmp_delay_idx)) *
          get_correct_factor_new_delay(
            curr_delays[kk], new_delays[kk], theta, group_idx, kk
          )['prob_proposing_new_value']
      }
    }
  }
  
  return(list(
    reverted_aug_dat = reverted_dat,
    correction_factor = correction_factor_rev
  ))
}

#' Calculate reverse proposal probability for E=0 to E=1 swaps
#' 
#' @description
#' Sequential step within `swap_Ei`. It performs the reverse of the
#' `perform_E0_to_E1_swap` step. It reverts the error indicators (`E`) and
#' dates (`D`) for the relevant dates back to their original state and
#' calculates the correction factor for this reverse move.
#' 
#' @details
#' Takes the partially reverted data from the previous step
#' (`reverse_resample_missing_dates`) and continues the reversal process. It
#' iterates backwards (`rev()`) through the dates that were originally `E=0`
#' and were swapped to `E=1` in the forward move.
#'
#' For each date, it reverts both the `E` and `D` values to their state in the
#' original data. It then calculates the log-probability of the forward move
#' from this newly reverted state back to the state before the revert. This
#' value is accumulated to create the total correction factor for this part of
#' the reverse path.
#' 
#' @param i Index of the individual.
#' @param group_idx Index of the group.
#' @param date_idx_E0_to_E1 Numeric vector of column indices for the dates
#'   that were `E=0` in the original data.
#' @param current_reverted_dat Augmented data list after the missing dates
#'   have already been reverted in the preceding step.
#' @param original_dat Augmented data list before any proposals were made.
#' @param theta List of model parameters (`mu`, `CV`, `zeta`).
#' @param obs_dat List of the observed data.
#' @param index_dates List defining the delays for each group.
#' @param range_dates Vector containing the overall range of dates.
#' 
#' @return A list containing two elements:
#' \itemize{
#'   \item{`reverted_aug_dat`}: Augmented data with the relevant dates
#'     reverted back to their original `E=0` state.
#'   \item{`correction_factor`}: Sum of the log-probabilities for this
#'     part of the reverse proposal.
#' }
#' @export
#' @seealso `perform_E0_to_E1_swap`, `compute_p_accept_move_E`
#' 
reverse_E0_to_E1_swap <- function(i,
                                  group_idx,
                                  date_idx_E0_to_E1,
                                  current_reverted_dat,
                                  original_dat,
                                  theta,
                                  obs_dat,
                                  index_dates,
                                  range_dates) {
  
  correction_factor_rev <- 0
  reverted_dat <- current_reverted_dat
  
  if (length(date_idx_E0_to_E1) > 0) {
    # Loop backwards
    for (k in rev(date_idx_E0_to_E1)) {
      state_before_revert <- reverted_dat
      
      # Revert E and D for this date to original values
      reverted_dat$E[[group_idx]][i, k] <- 0
      reverted_dat$D[[group_idx]][i, k] <- original_dat$D[[group_idx]][i, k]
      
      # Calculate reverse probability
      correction_factor_rev <- correction_factor_rev +
        compute_p_accept_move_E(direction = "E0_to_E1",
          i = i,
          group_idx = group_idx,
          date_idx = k,
          curr_aug_dat = state_before_revert,
          proposed_aug_dat = reverted_dat,
          theta = theta,
          obs_dat = obs_dat,
          index_dates = index_dates,
          range_dates = range_dates
        )[2]
    }
  }
  
  return(list(
    reverted_aug_dat = reverted_dat,
    correction_factor = correction_factor_rev
  ))
}

#' Identify delays affected by date changes
#' 
#' @param index_dates_group Index dates defining the delays for a group, where
#'  each column is a delay and the rows are the start/end date indices
#' @param all_changed_dates Unique numeric vector containing the indices of all
#'  dates that have been moved or changed e.g. from E0 to E1, from E1 to E0, and
#'  resampled if they were missing
#'  
#' @return Unique vector of the indices of all delays affected by the changed
#'  dates
#'
#' @export
#' 
#' @examples
#' # Example using 3 delays (1->2, 2->3, 1->4)
#' index_dates_group <- matrix(c(1, 2, 1,
#'                               2, 3, 4), nrow = 2, byrow = TRUE)
#' 
#' # 1. Change two dates in a hypothetical move
#' changed_dates <- c(2, 4)
#' 
#' # These dates are involved in all 3 delays
#' identify_affected_delays(index_dates_group, changed_dates)
#' 
#' # 2. Change one date in a hypothetical move
#' changed_dates <- 1
#' 
#' # Date 1 is only involved in delay indices 1 and 3
#' identify_affected_delays(index_dates_group, changed_dates)
#' 
#'
identify_affected_delays <- function(index_dates_group, all_changed_dates){
  
  # Return empty vector if no dates changed
  if (length(all_changed_dates) == 0) {
    return(integer(0))
  }
  
  # Find delays where any of the changed dates are involved
  affected_cols <- colSums(matrix(index_dates_group %in% all_changed_dates,
                                  nrow = nrow(index_dates_group))) > 0
  
  return(which(affected_cols))
}

#' Calculate the change in observation log-likelihood
#' 
#' @description
#' Calculate the change in the observation component of the log-likelihood
#' between a proposed and current state.
#' 
#' @param proposed_dat Proposed augmented data (e.g. `proposed_aug_dat_step3`)
#' @param current_dat Current augmented data (state before the new proposal)
#' @param i Index for individual
#' @param group_idx Index for group
#' @param all_changed_dates Vector of all date indices that were modified
#' @param obs_dat List of observed data
#' @param range_dates Vector of overall date range
#' 
#' @export
#' @return Numeric value of the log-ratio of observation likelihoods
#' 
calculate_obs_ratio <- function(proposed_dat, current_dat, i, group_idx,
                                all_changed_dates, obs_dat, range_dates) {
  
  ratio <- sum(
  LL_observation_term_by_group_delay_and_indiv(
    proposed_dat, obs_dat, group_idx,
    all_changed_dates, i, range_dates = range_dates
    ) - LL_observation_term_by_group_delay_and_indiv(
      current_dat, obs_dat, group_idx,
      all_changed_dates, i, range_dates = range_dates)
  )
  return(ratio)
}


#' Calculates the change in the error log-likelihood
#'
#' @description
#' Calculate the change in the error component of the log-likelihood
#' between a proposed and current state.
#' 
#' @param proposed_dat Proposed augmented data (e.g. `proposed_aug_dat_step3`)
#' @param current_dat Current augmented data (state before the new proposal)
#' @param i Index for individual
#' @param group_idx Index for group
#' @param all_changed_dates Vector of all date indices that were modified
#' @param theta List of model parameters, including `zeta`.
#' 
#' @export
#' @return Numeric value for the log-ratio of the error likelihoods
#' 
calculate_error_ratio <- function(proposed_dat, current_dat, i, group_idx,
                                  all_changed_dates, theta) {
  
  ratio <- sum(
    LL_error_term_by_group_delay_and_indiv(
      proposed_dat, theta, group_idx, all_changed_dates, i
    ) - LL_error_term_by_group_delay_and_indiv(
      current_dat, theta, group_idx, all_changed_dates, i
    )
  )
  return(ratio)
}

#' Calculates the change in the delay log-likelihood
#'
#' @description
#' Calculate the change in the delay component of the log-likelihood
#' between a proposed and current state.
#' 
#' @param proposed_dat Proposed augmented data (e.g. `proposed_aug_dat_step3`)
#' @param current_dat Current augmented data (state before the new proposal)
#' @param i Index for individual
#' @param group_idx Index for group
#' @param delay_idx Vector of delay indices affected by the date changes
#' @param theta List of model parameters, including `mu` and `CV`
#' @param obs_dat List of observed data
#' @param index_dates List defining the delays
#'
#' @export
#' @return Numeric value for the log-ratio of the delay likelihoods
#' 
calculate_delay_ratio <- function(proposed_dat, current_dat, i, group_idx,
                                  delay_idx, theta, obs_dat, index_dates) {
  
  ratio <- 0
  # Loop through each affected delay
  for (d in delay_idx) {
    ratio <- ratio + sum(
      LL_delays_term_by_group_delay_and_indiv(
        proposed_dat, theta, obs_dat, group_idx, d, i, index_dates
      ) - LL_delays_term_by_group_delay_and_indiv(
        current_dat, theta, obs_dat, group_idx, d, i, index_dates
      )
    )
  }
  return(ratio)
}


#' Performs one iteration of an MCMC move for the augmented data where the
#'  indicators of error in observations for one individual are swapped, i.e.
#'   the errors become non errors and vice versa.
#'
#' @param i The index of the individual for whom augmented data should be moved
#' @param group_idx The index of the group for whom augmented data should be
#'  moved. This should be an individual and for whom in curr_aug_dat there is at
#'   least one correctly recorded date and one erroneously recorded dates.
#' @param curr_aug_dat The current augmented data; a list of observed data, in
#'  the format returned by \code{\link{simul_true_data}}.
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#' @details \code{theta} should be a list containing:
#' \itemise{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups
#'   to be simulated data). Each element of \code{mu} should be a scalar of
#'    vector giving the mean delay(s) to use for simulation of dates in that
#'     group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar of vector giving the coefficient o variation
#'    of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
#' \itemise{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$CV}}
#' }
#' \code{index_dates} should be a list of length
#'  \code{n_groups = length(obs_dat)}. Each element of \code{index_dates}
#'   should be a matrix with 2 rows and a number of columns corresponding to
#'    the delays of interest for that group. For each column (i.e. each delay),
#'     the first row gives the index of the origin date, and the second row
#'      gives the index of the destination date.
#' The number of columns of index_dates[[k]] should match the length of
#'  theta$mu[[k]] and theta$CV[[k]]
#'
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and
#'  c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are
#'   respectively the mean and coefficient of variation of two delays: the
#'    first delay being between date 1 and date 2, and the second being between
#'     date 1 and date 3.
#'
#' The function performs the move as follows, using a Metropolis-Hastings
#'  algorithm.
#' It proposes to swap the indicators of errors in dates, e.g. 0 becomes 1 and
#'  1 becomes 0.
#' For the dates where we propose a move from E = 1 to E = 0, we automatically
#'  move D = the corresponding observed data.
#' For the dates where we propose a move from E = 0 to E = 1, D is then moved
#'  as follows: a new value is drawn from the marginal posterior of one of the
#'   delays this date is involved in, repeatedly until D falls on a different
#'    day than the observed date to be consistent with E = 1.
#' If the date is involved in several delays, one of the delays is randomly
#'  selected.
#' This move is not symmetrical so we use a correction factor in computing the
#'  probability of acceptance in the Metropolis Hastings which accounts for the
#'   asymmetry.
#' @return A list of two elements:
#'  \itemise{
#'  \item{\code{new_aug_dat}}{: Same as \code{curr_aug_dat} but where the
#'   relevant indicators of errors in dates have been updated}
#'  \item{\code{accept}}{: A scalar with value 1 if the move was accepted and 0
#'   otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
swap_Ei <- function(i,
                    group_idx,
                    curr_aug_dat,
                    theta,
                    obs_dat,
                    hyperparameters,
                    index_dates,
                    range_dates = NULL) {

  if (is.null(range_dates)) range_dates <- find_range(obs_dat)

  all_E_values <- curr_aug_dat$E[[group_idx]][i, ]
  date_idx <- seq_len(ncol(curr_aug_dat$E[[group_idx]]))

  date_idx_E0_to_E1 <- date_idx[all_E_values %in% 0]
  date_idx_E1_to_E0 <- date_idx[all_E_values %in% 1]
  date_idx_resample <- date_idx[all_E_values %in% -1]
  
  ## Forward moves ------------------------------------------------------------
  
  ## Step 1: Move E = 1 date(s) to E = 0
  proposed_aug_dat_step1 <- perform_E1_to_E0_swap(
    i, group_idx, date_idx_E1_to_E0,
    curr_aug_dat, theta, obs_dat,
    hyperparameters, index_dates, range_dates
    )
  
  ## Step 2: Move the original E = 0 dates to E = 1 ensuring that dates are
  ## plausible given the delay parameters
  step2_results <- perform_E0_to_E1_swap(
    i, group_idx, date_idx_E0_to_E1, proposed_aug_dat_step1,
    theta, obs_dat, hyperparameters, index_dates, range_dates
    )
  
  proposed_aug_dat_step2 <- step2_results$proposed_aug_dat
  corr_1_E0_to_E1 <- step2_results$correction_factor
  
  ## Step 3: Resample missing dates so that they are compatible with
  ## the new proposed dates.
  step3_results <- resample_missing_dates(
    i, group_idx, date_idx_resample, proposed_aug_dat_step2,
    theta, obs_dat, hyperparameters, index_dates, range_dates
  )
  
  proposed_aug_dat_step3 <- step3_results$proposed_aug_dat
  correct_factor_new_delay <- step3_results$correction_factor
  
  ## Reverse moves (for correction factor calc) -----------------------
  
  # Reverse of step 3
  rev_step3_results <- reverse_resample_missing_dates(
    i, group_idx, date_idx_resample,
    final_proposed_dat = proposed_aug_dat_step3,
    original_dat = curr_aug_dat, # the original state
    theta, index_dates
  )
  
  # Reverse of step 2
  rev_step2_results <- reverse_E0_to_E1_swap(
    i, group_idx, date_idx_E0_to_E1,
    current_reverted_dat = rev_step3_results$reverted_aug_dat,
    original_dat = curr_aug_dat,
    theta, obs_dat, index_dates, range_dates
  )

  # Identify all dates affected
  all_changed_dates <- unique(c(date_idx_E1_to_E0,
                                date_idx_E0_to_E1,
                                date_idx_resample))
  
  # Identify delays affected by date changes
  delay_idx <- identify_affected_delays(
    index_dates[[group_idx]],
    all_changed_dates
  )

  # Calculate the change in the observation likelihood term
  ratio_post_obs <- calculate_obs_ratio(proposed_dat = proposed_aug_dat_step3,
                                        current_dat = curr_aug_dat,
                                        i = i,
                                        group_idx = group_idx,
                                        all_changed_dates = all_changed_dates,
                                        obs_dat = obs_dat,
                                        range_dates = range_dates)

  # Calculate the change in the error likelihood term
  ratio_post_error <- calculate_error_ratio(proposed_dat = proposed_aug_dat_step3,
                                            current_dat = curr_aug_dat,
                                            i = i,
                                            group_idx = group_idx,
                                            all_changed_dates = all_changed_dates,
                                            theta = theta)

  # Calculate the change in the delay likelihood term
  ratio_post_delay <- calculate_delay_ratio(proposed_dat = proposed_aug_dat_step3,
                                            current_dat = curr_aug_dat,
                                            i = i,
                                            group_idx = group_idx,
                                            delay_idx = delay_idx,
                                            theta = theta,
                                            obs_dat = obs_dat,
                                            index_dates = index_dates)
  

  ratio_post <- ratio_post_obs + ratio_post_error + ratio_post_delay

  # This is not a symmetric move so need a correction factor

  # This correction factor should calculate
  # The probability of making the 3 types of move (E1 to E0, E0 to E1, and NAs) sequentially
  # and the probability of making the reverse move (in that same order)

  ### First type of move
  # Each move from E1 to E0 is made with probability 1 so that's easy (0 correction on the log scale)

  ### Second type of move
  # Each move from E0 to E1 is made with probability

  ### Third type of move
  # Moving missing data
  # Calculation of correction factors is made above

  ### Recap
  # probability of moves from E1 to E0 is 1 so corresponding term is zero
  # corr_1_E0_to_E1 # computes (- log_p_proposing_move_E0old_to_E1new)
  # corr_1_E0_to_E1_rev # computes (- log_p_proposing_move_from_E0new_to_E1old)
  # correct_factor_new_delay # computes the log_prob_proposing_new_missing_values
  # correct_factor_new_delay_rev # computes the log_prob_proposing_old_missing_values

  corr <- (rev_step2_results$correction_factor - step2_results$correction_factor) +
          (rev_step3_results$correction_factor - step3_results$correction_factor)

  #print(paste(c(ratio_post ,
  #                corr_1_E0_to_E1 , corr_1_E0_to_E1_rev ,
  #              correct_factor_new_delay_rev , correct_factor_new_delay)))

  # if(is.infinite(ratio_post) & is.infinite(corr)) {
  #   browser()
  #   # this seems to yield problems currently when two dates are wrong and one right
  # }

  p_accept <- ratio_post + corr

  res <- decide_acceptance(log_p_accept = p_accept,
                           proposed_dat = proposed_aug_dat_step3,
                           current_dat = curr_aug_dat)

  return(res)

}

# ----------------------------------------------------------------------------
# Move mean or CV of delay
# ----------------------------------------------------------------------------

#' Performs one iteration of an MCMC move for either the parameter mu or the
#'  parameter CV (mean or CV of the various delays to be estimated)
#'
#' @param what A string ("mu" or "CV") indicating which of the parameters to
#'  move
#' @param group_idx The index of the group for which mu or CV should be moved
#' @param delay_idx The index of the delay for which mu or CV should be moved
#' @param sdlog The standard deviation to be used for the proposal distribution,
#'  see details.
#' @param aug_dat The augmented data; a list of observed data, in the format
#'  returned by \code{\link{simul_true_data}}.
#' @param curr_theta The current list of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @details \code{curr_theta} should be a list containing:
#' \itemise{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups
#'   to be simulated data). Each element of \code{mu} should be a scalar of
#'    vector giving the mean delay(s) to use for simulation of dates in that
#'     group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar of vector giving the coefficient o variation
#'    of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
#' \itemise{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape
#'   parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$CV}}
#' }
#'
#' \code{index_dates} should be a list of length
#'  \code{n_groups = length(obs_dat)}. Each element of \code{index_dates}
#'   should be a matrix with 2 rows and a number of columns corresponding to
#'    the delays of interest for that group. For each column (i.e. each delay),
#'     the first row gives the index of the origin date, and the second row
#'      gives the index of the destination date.
#' The number of columns of index_dates[[k]] should match the length of
#'  theta$mu[[k]] and theta$CV[[k]]
#'
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and
#'  c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are
#'   respectively the mean and coefficient of variation of two delays: the
#'    first delay being between date 1 and date 2, and the second being between
#'     date 1 and date 3.
#'
#' The function performs the move as follows, using a Metropolis-Hastings
#'  algorithm.
#' For the parameter to be moved, a new value is drawn from a lognormal
#'  distribution with parameters \code{meanlog} equals the log of the current
#'   parameter value, and \code{sdlog=sdlog}.
#' The new parameter set is then accepted with probability given by the ratio
#'  of the posterior values at the new parameter set and the old parameter set,
#'   multiply by a correction factor to reflect the non-symetrical nature of
#'    the move.
#' @return A list of two elements:
#'  \itemise{
#'  \item{\code{new_theta}: Same as \code{curr_theta} but where
#'   \code{curr_theta$zeta} has been updated}
#'  \item{\code{accept}: A scalar with value 1 if the move was accepted and 0
#'   otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
move_lognormal <- function(what = c("mu", "CV"),
                           group_idx,
                           delay_idx,
                           sdlog,
                           aug_dat,
                           curr_theta,
                           obs_dat,
                           hyperparameters,
                           index_dates) {

  what <- match.arg(what)

  # draw proposed value
  curr_param_value <- curr_theta[[what]][[group_idx]][delay_idx]
  proposed_param_value <- rlnorm(
    1, meanlog = log(curr_param_value), sdlog = sdlog
  )

  proposed_theta <- curr_theta
  proposed_theta[[what]][[group_idx]][delay_idx] <- proposed_param_value

  # calculates probability of acceptance
  if (what == "mu") {
    ratio_post <- lprior_params_delay(what, proposed_theta, hyperparameters) -
      lprior_params_delay(what, curr_theta, hyperparameters)
  } else if (what == "CV") {
    ratio_post <- lprior_params_delay(what, proposed_theta, hyperparameters) -
      lprior_params_delay(what, curr_theta, hyperparameters)
  }

  delta <- compute_delta_group_delay_and_indiv(
    aug_dat$D, group_idx, seq_len(nrow(obs_dat[[group_idx]])),
    delay_idx, index_dates
  ) # same for proposed and current par values so no need to recompute twice

  ratio_post <- ratio_post + sum(
    LL_delays_term_by_group_delay_and_indiv(
      aug_dat, proposed_theta, obs_dat, group_idx, delay_idx,
      seq_len(nrow(obs_dat[[group_idx]])), index_dates, delta
    )
  ) - sum(
    LL_delays_term_by_group_delay_and_indiv(
      aug_dat, curr_theta, obs_dat, group_idx, delay_idx,
      seq_len(nrow(obs_dat[[group_idx]])), index_dates, delta
    )
  )

  ### note that ratio_post should be the same as:
  # ratio_post_long <- lposterior_total(aug_dat, proposed_theta, obs_dat,
  # hyperparameters, index_dates) -
  # lposterior_total(aug_dat, curr_theta, obs_dat, hyperparameters, index_dates)

  # correction for lognormal distribution
  correction <- log(proposed_param_value) - log(curr_param_value)
  # things are additive here as on log scale:
  p_accept <- ratio_post + correction
  if (p_accept > 0) p_accept <- 0

  # accept/reject step
  tmp <- log(runif(1))
  if (tmp < p_accept) { # accepting with a certain probability
    new_theta <- proposed_theta
    accept <- 1
  } else { # reject
    new_theta <- curr_theta
    accept <- 0
  }

  return(list(new_theta = new_theta, accept = accept))

}
# test_move_mu <- move_lognormal(what="mu", group_idx=1, delay_idx=1, sdlog=0.1,
# aug_dat, curr_theta = theta, obs_dat, hyperparameters)
# test_move_mu$new_theta$mu[[1]][1] # new value
# theta$mu[[1]][1] # old value


# ----------------------------------------------------------------------------
# Move zeta (probability of erroneous recording of dates)
# ----------------------------------------------------------------------------

#' Performs one iteration of an MCMC move for the parameter zeta (probability
#'  of a data being recorded erroneously, given it is recorded)
#'
#' @param aug_dat The augmented data; a list of observed data, in the format
#'  returned by \code{\link{simul_true_data}}.
#' @param curr_theta The current list of parameters; see details.
#' @param hyperparameters A list of hyperparameters: see details.
#' @details \code{curr_theta} should be a list containing:
#' \itemise{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups
#'   to be simulated data). Each element of \code{mu} should be a scalar of
#'    vector giving the mean delay(s) to use for simulation of dates in that
#'     group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar of vector giving the coefficient o variation
#'    of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
#' \itemise{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape
#'   parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$CV}}
#' }
#'
#' The function performs the move, using a Gibbs sampler.
#' A new value of parameter zeta is drawn from its marginal posterior
#'  distribution, that is a beta distribution with parameters:
#' \itemise{
#'  \item{\code{first shape parameter}}{: Equal to
#'   \code{hyperparameters$shape1_prob_error} + number_of_errors, where
#'    number_of_errors is the number of data points recorded with errors}
#'  \item{\code{second shape parameter}}{: Equal to
#'   \code{hyperparameters$shape2_prob_error} + number_of_recorded_dates -
#'    number_of_errors, where number_of_recorded_dates is the number of data
#'     points which are not missing, and number_of_errors is the number of data
#'      points recorded with errors}
#' }
#' @return A list of two elements:
#'  \itemise{
#'  \item{\code{new_theta}: Same as \code{curr_theta} but where
#'   \code{curr_theta$zeta} has been updated}
#'  \item{\code{accept}: A scalar with value 1 (as we are using a Gibbs
#'   sampler the move is always accepted)}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
move_zeta_gibbs <- function(aug_dat,
                            curr_theta,
                            hyperparameters) {

  tmp <- compute_n_errors(aug_dat)
  number_of_errors <- tmp[1]
  number_of_recorded_dates <- tmp[2]

  # drawing from the marginal posterior distribution directly
  new_zeta <- rbeta(
    1,
    shape1 = hyperparameters$shape1_prob_error + number_of_errors,
    shape2 = hyperparameters$shape2_prob_error +
      number_of_recorded_dates - number_of_errors
  )

  # therefore accept automatically
  new_theta <- curr_theta
  new_theta$zeta <- new_zeta
  accept <- 1

  list(new_theta = new_theta, accept = accept)

}
# test_move_zeta_gibbs <- move_zeta_gibbs(aug_dat, theta, hyperparameters)
# test_move_zeta_gibbs$new_theta$zeta # new value
# theta$zeta # old value
