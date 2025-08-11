###############################################
# Move functions
###############################################

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
#' \itemize{
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
#' \itemize{
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
#'  \itemize{
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
    proposed_aug_dat, theta, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  )

  ll_current <- LL_observation_term_by_group_delay_and_indiv(
    curr_aug_dat, theta, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  )

  ratio_post <- ll_proposed - ll_current

  # Add error term difference only if E changed
  different_E <- proposed_aug_dat$E[[group_idx]][i, date_idx] !=
    curr_aug_dat$E[[group_idx]][i, date_idx]

  if (any(different_E)) {
    ratio_post <- ratio_post +
      LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat,
                                             group_idx, date_idx, i) -
      LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat,
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
  if (p_accept > 0) p_accept <- 0

  # Accept/reject the proposal ------------------------------------------------
  if (log(runif(1)) < p_accept) {
    new_aug_dat <- proposed_aug_dat
    accept <- 1
  } else {
    new_aug_dat <- curr_aug_dat
    accept <- 0
  }

  return(list(new_aug_dat = new_aug_dat, accept = accept))

}
# test_move_Di <- move_Di(i=1, group_idx=1, date_idx=1, curr_aug_dat = aug_dat,
# theta, obs_dat, hyperparameters)
# test_move_Di$new_aug_dat$D[[1]][1,1] # new value
# aug_dat$D[[1]][1,1] # old value


# -----------------------------------------------------------------------------
# Move augmented indicator for whether date is correctly recorded E
# -----------------------------------------------------------------------------

#' Propose a new value for the true date when transitioning the error indicator
#'  from 0 (no error) to 1 (observed with error), ensuring the new date is
#'   different from the observed date.
#'
#' @param i Index of individual(s) for whom augmented data should be moved.
#' @param group_idx Index of the group for whom augmented data should be moved.
#' @param date_idx Index of the date which should be moved.
#' @param curr_aug_dat The current augmented data.
#' @param theta List of parameters, including mu, CV and zeta.
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters.
#' @param index_dates A list containing the delays defined for each group.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#'
#' @return A vector of proposed true dates that differ from the observed dates.

## ANNE: TODO: clarify in comments that this actually doesn't move E just the corresponding D.
propose_move_from_E0_to_E1 <- function(i,
                                       group_idx,
                                       date_idx,
                                       curr_aug_dat,
                                       theta,
                                       obs_dat,
                                       hyperparameters,
                                       index_dates,
                                       range_dates = NULL) {

  # ANNE: which delays is this date / are these dates involved in
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

  # ANNE: retrieve the corresponding delay parameters
  param_delay <- lapply(seq_along(date_idx), function(e) {
    find_params_gamma(
      theta$mu[[group_idx]][which_delay[[e]]],
      CV = theta$CV[[group_idx]][which_delay[[e]]]
    )
  })

  # ANNE: store the current corresponding augmented dates
  # ANNE: TODO: remove this line, not needed in this function
  curr_aug_dat_value <- curr_aug_dat$D[[group_idx]][i, date_idx]

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

#' Computes acceptance probability for a move from E = 0 (date observed
#'  correctly) to E = 1 (date observed incorrectly).
#'
#' @param i Index of individual(s) for whom augmented data should be moved.
#' @param group_idx Index of the group for whom augmented data should be moved.
#' @param date_idx Index of the date which should be moved.
#' @param curr_aug_dat The current augmented data.
#' @param proposed_aug_dat Proposed augmented data.
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#'
#'  @return Vector of length 2. The first element is the difference in log
#'   posterior between the proposed and current augmented data. The second
#'    is the log proposal correction factor needed to adjust for asymmetry in
#'     the proposal distribution.
#'
compute_p_accept_move_from_E0_to_E1 <- function(i,
                                                group_idx,
                                                date_idx,
                                                curr_aug_dat,
                                                proposed_aug_dat,
                                                theta,
                                                obs_dat,
                                                hyperparameters,
                                                index_dates,
                                                range_dates = NULL) {

  proposed_aug_dat_value <- proposed_aug_dat$D[[group_idx]][i, date_idx]

  # Delays that are affected by the change in date date_idx
  delay_idx <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)[, 2]

  # Compute log posterior difference (proposed - current) ---------------------

  # Difference in observation likelihood
  ratio_post <- LL_observation_term_by_group_delay_and_indiv(
    proposed_aug_dat, theta, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  ) - LL_observation_term_by_group_delay_and_indiv(
    curr_aug_dat, theta, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  )

  # Difference in error likelihood
  ratio_post <- ratio_post + LL_error_term_by_group_delay_and_indiv(
    proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i
  ) - LL_error_term_by_group_delay_and_indiv(
    curr_aug_dat, theta, obs_dat, group_idx, date_idx, i
  )

  # For each affected delay, difference in delay likelihood
  for (d in delay_idx) {
    ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(
      proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates
    ) - LL_delays_term_by_group_delay_and_indiv(
      curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates
    )
  }

  # Combine all log likelihood differences
  ratio_post <- sum(ratio_post)

  ### note that ratio_post should be the same as:
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat,
  # hyperparameters, index_dates) -
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates)

  # ANNE: the above now seems to work but maybe worth checking a bit more thoroughly in tests

  # Correct asymmetry ---------------------------------------------------------

  # Index for the date within index_dates
  x <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)

  # Index for the delay
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

  # Proposal correction factor for each delay
  find_correction_factor <- function(e) {
    if (date_idx < from_idx[e]) {
      delay <- from_value[e] - proposed_aug_dat_value
      forbidden_delay <- from_value[e] - obs_dat[[group_idx]][i, date_idx] # ANNE: this is the delay corresponding to the observed date, hence would keep E = 0 and not allow a move to E = 1
    } else {
      delay <- proposed_aug_dat_value - from_value[e]
      forbidden_delay <- obs_dat[[group_idx]][i, date_idx] - from_value[e]
    }

    # Mean and CV of delays
    mu <- theta$mu[[group_idx]][which_delay[e]]
    cv <- theta$CV[[group_idx]][which_delay[e]]

    # Probability mass for delay after adjusting for invalid delay
    K <- DiscrGamma(k = delay, mu = mu, cv = cv, log = FALSE) / # ANNE: calculating the probability of randomly drawing the proposed delay
      (1 - DiscrGamma(k = forbidden_delay, mu = mu, cv = cv, log = FALSE)) # ANNE: renormalising because we do not allow this specific delay

    K # probability of drawing the proposed delay if this specific delay index is chosen
  }

  # Average correction factors over all delays involving this date
  # ANNE: this accounts for the fact that the proposal distribution is a mixture distribution
  # where the index of the delay used for sampling is drawn at random.
  # hence the correction factor, which should be the probability of drawing the proposed delay irrespective of which delay index was chosen
  # hence it's P(proposed delay | delay index 1) * P(delay index 1) + P(proposed delay | delay index 2) * P(delay index 2) + ...
  # which simplifies to just the mean of P(proposed delay | delay index i) because all the P(delay index i) are the same
  K <- mean(sapply(seq_along(which_delay), find_correction_factor))

  # Log proposal correction factor
  # ANNE: this should be calculated as logP(proposing current delay) - logP(proposing new delay)
  # but this: logP(proposing current delay) is zero because P(proposing current delay) = 1 because there is only 1 way of moving to the observed date i.e. to generate E = 0
  # hence logcorrection <- 0 - logP(proposing new delay)
  logcorrection <- -log(K) # log_p_move_from_new_to_old - log_p_move_from_old_to_new

  return(c(ratio_post, logcorrection))
}


#' Propose the observed date when transitioning the error indicator from 1
#'  (observed with error) to 0 (observed with no error). This move assumes the
#'  true date was correctly recorded and simply sets D to the observed date.
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
#' @return Proposed true date, which is the same as the observed date.

propose_move_from_E1_to_E0 <- function(i,
                                       group_idx,
                                       date_idx,
                                       curr_aug_dat, # remove
                                       theta, # remove
                                       obs_dat,
                                       hyperparameters, # remove
                                       index_dates, # remove
                                       range_dates) { # remove

  proposed_aug_dat_value <- obs_dat[[group_idx]][i, date_idx]

  proposed_aug_dat_value
}


#' Calculates the log acceptance probability for an MCMC move that changes the
#'  error indicator from E = 1 (observed with error) to E = 0 (observed without
#'   error) and sets the true date to the observed date.
#'
#' @param i Index of individual(s) for whom augmented data should be moved.
#' @param group_idx Index of the group for whom augmented data should be moved.
#' @param date_idx Index of the date which should be moved.
#' @param curr_aug_dat The current augmented data.
#' @param proposed_aug_dat The proposed augmented data.
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#'
#'  @return A vector of length 2. The first element is the log posterior
#'   difference (proposed - current), and the second element is the proposal
#'    asymmetry correction.

compute_p_accept_move_from_E1_to_E0 <- function(i,
                                                group_idx,
                                                date_idx,
                                                curr_aug_dat,
                                                proposed_aug_dat,
                                                theta,
                                                obs_dat,
                                                hyperparameters,
                                                index_dates,
                                                range_dates) {

  ### ANNE: this is the exact opposite move from compute_p_accept_move_from_E0_to_E1
  ## hence most of the code is the same
  ## the ratio of the posteriors is just as usual Post(new proposed value) - Post (old value)
  ## and the probability of accepting a move is calculated exactly in the opposite way compared to compute_p_accept_move_from_E0_to_E1
  ## i.e. there is only 1 way of moving from E1 to E0 (because there is only one date that is equal to the observed date), so P(proposing the new value) = 1 and hence the log is zero
  ## and the Probability of moving from E0 to this specific E1 and D combination is calculated according to the delay distribution, i.e. DircGamma, but discounting the one value we cannot choose
  ## because it corresponds to E = 0 not E = 1.

  # Current date before the move
  curr_aug_dat_value <- curr_aug_dat$D[[group_idx]][i, date_idx]

  # Delays that are affected by the change in date date_idx
  delay_idx <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)[, 2]

  # Compute log posterior difference (proposed - current) ---------------------

  # Difference in observation likelihood
  ratio_post <- LL_observation_term_by_group_delay_and_indiv(
    proposed_aug_dat, theta, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  ) - LL_observation_term_by_group_delay_and_indiv(
    curr_aug_dat, theta, obs_dat,
    group_idx, date_idx, i, range_dates = range_dates
  )

  # Difference in error likelihood
  ratio_post <- ratio_post + LL_error_term_by_group_delay_and_indiv(
    proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i
  ) - LL_error_term_by_group_delay_and_indiv(
    curr_aug_dat, theta, obs_dat, group_idx, date_idx, i
  )

  # For each affected delay, difference in delay likelihood
  for (d in delay_idx)
    ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(
      proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates
    ) - LL_delays_term_by_group_delay_and_indiv(
      curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates
    )

  # Combine all log-likelihood differences
  ratio_post <- sum(ratio_post)

  ### note that ratio_post should be the same as:
  ## ANNE: TODO check this works
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat,
  # hyperparameters, index_dates) -
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates)

  # Correct asymmetry ---------------------------------------------------------

  # Index for the date within index_dates
  x <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)

  # Index for the delay
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

  # Proposal correction factor for each delay
  find_correction_factor_2 <- function(e) {
    ## ANNE: TODO check if this is any different from the find_correction_factor function,
    ## if so perhaps move out of these move functions so you can merge into one single function
    if (date_idx < from_idx[e]) {
      delay <- from_value[e] - curr_aug_dat_value
      forbidden_delay <- from_value[e] - obs_dat[[group_idx]][i, date_idx]
    } else {
      delay <- curr_aug_dat_value - from_value[e]
      forbidden_delay <- obs_dat[[group_idx]][i, date_idx] - from_value[e]
    }

    # Mean and CV of delays
    mu <- theta$mu[[group_idx]][which_delay[e]]
    cv <- theta$CV[[group_idx]][which_delay[e]]

    # Probability mass for delay adjusting for invalid delay
    K <- DiscrGamma(k = delay, mu = mu, cv = cv, log = FALSE) /
      (1 - DiscrGamma(k = forbidden_delay, mu = mu, cv = cv, log = FALSE))

    K
  }

  # Average correction factors over all delays involving this date
  K <- mean(sapply(seq_along(which_delay), find_correction_factor_2))

  # Log proposal correction factor
  logcorrection <- +log(K)

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
#' \itemize{
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
#' \itemize{
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
#'  \itemize{
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
    warning(
      "In move_Ei, i should be a numeric, not a vector. Using i[1] instead."
    )
  }

  if (is.null(range_dates)) range_dates <- find_range(obs_dat)

  curr_E_value <- curr_aug_dat$E[[group_idx]][i, date_idx]
  proposed_aug_dat <- curr_aug_dat

  # if data not missing
  if (curr_E_value != -1) {

    new_E_value <- 1 - curr_E_value
    proposed_aug_dat$E[[group_idx]][i, date_idx] <- new_E_value

    # moving from E=0 to E=1
    if (curr_E_value == 0) {

      proposed_aug_dat$D[[group_idx]][i, date_idx] <- propose_move_from_E0_to_E1(
        i, group_idx, date_idx, curr_aug_dat, theta,
        obs_dat, hyperparameters, index_dates, range_dates
      )

      tmp <- compute_p_accept_move_from_E0_to_E1(
        i, group_idx, date_idx, curr_aug_dat, proposed_aug_dat,
        theta, obs_dat, hyperparameters, index_dates, range_dates
      )

      if (any(is.infinite(tmp))) p_accept <- -Inf else p_accept <- sum(tmp)
      if (p_accept > 0) p_accept <- 0

      # accept/reject step
      tmp <- log(runif(1))
      if (tmp < p_accept) { # accepting with a certain probability
        new_aug_dat <- proposed_aug_dat
        accept <- 1
      } else { # reject
        new_aug_dat <- curr_aug_dat
        accept <- 0
      }

      # return a list of size 2 where
      #		the first value is the new augmented data set in the chain
      #		the second value is 1 if the proposed value was accepted, 0 otherwise
      res <- list(new_aug_dat = new_aug_dat, accept = accept)

    } else if (curr_E_value == 1) { # moving from E=1 to E=0

      proposed_aug_dat$D[[group_idx]][i, date_idx] <- propose_move_from_E1_to_E0(
        i, group_idx, date_idx, curr_aug_dat, theta,
        obs_dat, hyperparameters, index_dates, range_dates
      )

      tmp <- compute_p_accept_move_from_E1_to_E0(
        i, group_idx, date_idx, curr_aug_dat, proposed_aug_dat,
        theta, obs_dat, hyperparameters, index_dates, range_dates
      )

      if (any(is.infinite(tmp))) p_accept <- -Inf else p_accept <- sum(tmp)
      if (p_accept > 0) p_accept <- 0

      # accept/reject step
      tmp <- log(runif(1))
      if (tmp < p_accept) { # accepting with a certain probability
        new_aug_dat <- proposed_aug_dat
        accept <- 1
      } else { # reject
        new_aug_dat <- curr_aug_dat
        accept <- 0
      }

      # return a list of size 2 where
      #		the first value is the new augmented data set in the chain
      #		the second value is 1 if the proposed value was accepted, 0 otherwise
      res <- list(new_aug_dat = new_aug_dat, accept = accept)
    }
  } else { # if E=-1, can't move.
    res <- list(curr_aug_dat = curr_aug_dat, accept = 0)
  }
  return(res)
}

# ----------------------------------------------------------------------------
# Swap Es - when only 2 Es, related by a delay, are recorded, one with error
# and one without error, propose to swap the two
# ----------------------------------------------------------------------------

# Identify rows where (other than -1s) errors are not all the same
# i.e. mixed errors and non-errors
## ANNE: TODO for Rebecca to add documentation for this function
## ANNE: this function finds individuals that have among all their dates one that is E = 0 and one that is E = 1 which hence could be swapped
find_Eis_to_swap <- function(group_idx, curr_aug_dat) {
  Es <- curr_aug_dat$E[[group_idx]]
  # check for more than one unique non_missing entry
  which(sapply(seq_len(nrow(Es)), function(i) {
    non_missing <- Es[i, Es[i, ] != -1]
    unique_values_of_E_i <- unique(non_missing)
    length(unique_values_of_E_i) > 1 # this only has length >1 if it contains both a zero and a one.
  }))
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
#' \itemize{
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
#' \itemize{
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
#'  \itemize{
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

  # TO DO: Add tests to see how these handle multiple values
  date_idx_E0_to_E1 <- date_idx[all_E_values %in% 0]
  date_idx_E1_to_E0 <- date_idx[all_E_values %in% 1]
  date_idx_resample <- date_idx[all_E_values %in% -1]
  
  ## Step 1: Move E = 1 date(s) to E = 0
  proposed_aug_dat_step1 <- curr_aug_dat
  if (length(date_idx_E1_to_E0) > 0) {
    for (k in date_idx_E1_to_E0) {
      proposed_aug_dat_step1$E[[group_idx]][i, k] <- 0
      proposed_aug_dat_step1$D[[group_idx]][i, k] <-
        propose_move_from_E1_to_E0(
          i, group_idx, k, curr_aug_dat, theta, obs_dat,
          hyperparameters, index_dates, range_dates
        )
    }
  }
  
  ## Step 2: Move the original E = 0 dates to E = 1 ensuring that dates are
  ## plausible given the delay parameters
  proposed_aug_dat_step2 <- proposed_aug_dat_step1
  corr_1_E0_to_E1 <- 0
  if (length(date_idx_E0_to_E1) > 0) {
    for (k in date_idx_E0_to_E1) {
      current_state <- proposed_aug_dat_step2
      proposed_aug_dat_step2$E[[group_idx]][i, k] <- 1
      proposed_aug_dat_step2$D[[group_idx]][i, k] <-
        propose_move_from_E0_to_E1(
          i, group_idx, k, current_state,
          theta, obs_dat, hyperparameters, index_dates, range_dates
        )
      corr_1_E0_to_E1 <- corr_1_E0_to_E1 +
        compute_p_accept_move_from_E0_to_E1(
          i = i,
          group_idx = group_idx,
          date_idx = k,
          curr_aug_dat = current_state,
          proposed_aug_dat = proposed_aug_dat_step2,
          theta = theta,
          obs_dat = obs_dat,
          hyperparameters = hyperparameters,
          index_dates = index_dates,
          range_dates = range_dates
        )[2]
    }
  }

  ## Step 3: Resample missing dates so that they are compatible with
  ## the new proposed dates.
  proposed_aug_dat_step3 <- proposed_aug_dat_step2
  correct_factor_new_delay <- 0
  if (length(date_idx_resample) > 0) {
    for (k in date_idx_resample) {
      tmp_delay <- propose_new_delay(
        i, group_idx, k, proposed_aug_dat_step3,
        theta, obs_dat, hyperparameters, index_dates, range_dates
      )
      # propose a new date compatible with newly drawn dates
      proposed_aug_dat_step3 <- tmp_delay$proposed_aug_dat
      correct_factor_new_delay <- correct_factor_new_delay +
        get_correct_factor_new_delay(
          tmp_delay$curr_delay, tmp_delay$sample_delay, theta,
          group_idx, tmp_delay$which_delay
        )["prob_proposing_new_value"]
    }
  }

  ## Mimic the reverse moves for correction factor calc
  # # step 1
  # proposed_aug_dat_rev1 <- proposed_aug_dat_step3
  # proposed_aug_dat_rev1$E[[group_idx]][i, date_idx_E0_to_E1] <- 0
  # proposed_aug_dat_rev1$D[[group_idx]][i, date_idx_E0_to_E1] <-
  #   propose_move_from_E1_to_E0(i, group_idx, date_idx_E0_to_E1,
  #                              proposed_aug_dat_step3, theta, obs_dat, hyperparameters,
  #                              index_dates, range_dates)
  
  
  # Reverse of step 3
  correct_factor_new_delay_rev <- 0
  proposed_aug_dat_rev3 <- proposed_aug_dat_step3
  curr_delays <- compute_delta(proposed_aug_dat_rev3$D, index_dates)[[group_idx]][i, ]
  if(length(date_idx_resample) > 0) {
    for(k in date_idx_resample) {
      proposed_aug_dat_rev3$D[[group_idx]][i, k] <- curr_aug_dat$D[[group_idx]][i, k]
    }
    new_delays <- compute_delta(proposed_aug_dat_rev3$D, index_dates)[[group_idx]][i, ]
    
    for(k in date_idx_resample) {
      tmp_delay_idx <- which(index_dates[[group_idx]] == k, arr.ind = TRUE)[, 2]
      for(kk in tmp_delay_idx) {
        correct_factor_new_delay_rev <- correct_factor_new_delay_rev +
          (1/length(tmp_delay_idx)) *
          get_correct_factor_new_delay(
            curr_delays[kk], new_delays[kk], theta, group_idx, kk
          )['prob_proposing_new_value']
      }
    }
  }
  
  # Reverse of step 2
  corr_1_E0_to_E1_rev <- 0
  proposed_aug_dat_rev2 <- proposed_aug_dat_rev3
  
  if (length(date_idx_E0_to_E1) > 0) {
    for (k in rev(date_idx_E0_to_E1)) {
      current_state <- proposed_aug_dat_rev2
      proposed_aug_dat_rev2$E[[group_idx]][i, k] <- 0
      proposed_aug_dat_rev2$D[[group_idx]][i, k] <- curr_aug_dat$D[[group_idx]][i, k]
      
      corr_1_E0_to_E1_rev <- corr_1_E0_to_E1_rev +
        compute_p_accept_move_from_E0_to_E1(
          i = i,
          group_idx = group_idx,
          date_idx = k,
          curr_aug_dat = current_state,
          proposed_aug_dat = proposed_aug_dat_rev2,
          theta = theta,
          obs_dat = obs_dat,
          hyperparameters = hyperparameters,
          index_dates = index_dates,
          range_dates = range_dates
        )[2]
    }
  }
  

  delay_idx <- which(
    colSums(matrix(index_dates[[group_idx]] %in% date_idx_E1_to_E0,
                   nrow = nrow(index_dates[[group_idx]]))) > 0
  )

  delay_idx <- c(delay_idx, which(
    colSums(matrix(index_dates[[group_idx]] %in% date_idx_E0_to_E1,
                   nrow = nrow(index_dates[[group_idx]]))) > 0
  ))

  if(length(date_idx_resample) > 0) {
    delay_idx <- c(delay_idx, which(
      colSums(matrix(index_dates[[group_idx]] %in% date_idx_resample,
                     nrow = nrow(index_dates[[group_idx]]))) > 0
    ))
  }

  delay_idx <- sort(unique(delay_idx)) ## ANNE: TODO: I think by definition this will be all the delays!

  ratio_post_obs <- sum(
    LL_observation_term_by_group_delay_and_indiv(
      proposed_aug_dat_step3, theta, obs_dat, group_idx,
      date_idx_E1_to_E0, i, range_dates = range_dates
    ) - LL_observation_term_by_group_delay_and_indiv(
      curr_aug_dat, theta, obs_dat, group_idx,
      date_idx_E1_to_E0, i, range_dates = range_dates)
  ) + sum(
    LL_observation_term_by_group_delay_and_indiv(
      proposed_aug_dat_step3, theta, obs_dat, group_idx,
      date_idx_E0_to_E1, i, range_dates = range_dates
    ) - LL_observation_term_by_group_delay_and_indiv(
      curr_aug_dat, theta, obs_dat, group_idx,
      date_idx_E0_to_E1, i, range_dates = range_dates)
  )

  if(length(date_idx_resample) > 0) {
    ratio_post_obs <- ratio_post_obs + sum(
      LL_observation_term_by_group_delay_and_indiv(
        proposed_aug_dat_step3, theta, obs_dat, group_idx,
        date_idx_resample, i, range_dates = range_dates
      ) - LL_observation_term_by_group_delay_and_indiv(
        curr_aug_dat, theta, obs_dat, group_idx,
        date_idx_resample, i, range_dates = range_dates)
    )
  }
  ## should be the same as:
  # LL_observation_term(proposed_aug_dat_step3, theta, obs_dat, range_dates) -
  # LL_observation_term(curr_aug_dat, theta, obs_dat, range_dates)

  ratio_post_error <- sum(
    LL_error_term_by_group_delay_and_indiv(
      proposed_aug_dat_step3, theta, obs_dat, group_idx, date_idx_E1_to_E0, i
    ) - LL_error_term_by_group_delay_and_indiv(
      curr_aug_dat, theta, obs_dat, group_idx, date_idx_E1_to_E0, i
    )
  ) + sum(
    LL_error_term_by_group_delay_and_indiv(
      proposed_aug_dat_step3, theta, obs_dat, group_idx, date_idx_E0_to_E1, i
    ) - LL_error_term_by_group_delay_and_indiv(
      curr_aug_dat, theta, obs_dat, group_idx, date_idx_E0_to_E1, i
    )
  )
  if(length(date_idx_resample) > 0) {
    ratio_post_error <- ratio_post_error + sum(
      LL_error_term_by_group_delay_and_indiv(
        proposed_aug_dat_step3, theta, obs_dat, group_idx, date_idx_resample, i
      ) - LL_error_term_by_group_delay_and_indiv(
        curr_aug_dat, theta, obs_dat, group_idx, date_idx_resample, i
      )
    )
  }
  ## should be the same as:
  # LL_error_term(proposed_aug_dat_step3, theta, obs_dat) -
  # LL_error_term(curr_aug_dat, theta, obs_dat)
  ## LL_error_term_slow(proposed_aug_dat_step3, theta, obs_dat) -
  ## LL_error_term_slow(curr_aug_dat, theta, obs_dat)

  ratio_post_delay <- 0
  for (d in delay_idx) {
    ratio_post_delay <- ratio_post_delay + sum(
      LL_delays_term_by_group_delay_and_indiv(
        proposed_aug_dat_step3, theta, obs_dat, group_idx, d, i, index_dates
      ) - LL_delays_term_by_group_delay_and_indiv(
        curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates
      )
    )
  }
  ## should be the same as:
  # LL_delays_term(proposed_aug_dat_step3, theta, obs_dat, index_dates) -
  # LL_delays_term(curr_aug_dat, theta, obs_dat, index_dates)

  ratio_post <- ratio_post_obs + ratio_post_error + ratio_post_delay

  ### should be the same as:
  ## ANNE: this may need to be checked more thoroughly but does work on one example
  # ratio_post_long <- lposterior_total(proposed_aug_dat_step3, theta, obs_dat,
  # hyperparameters, index_dates, range_dates) -
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates, range_dates)

  # This is not a symmetric move so need a correction factor

  # This correction factor should calculate
  # The probability of making the 3 types of move (E1 to E0, E0 to E1, and NAs) sequentially
  # and the probability of making the reverse move (in that same order)

  ### First type of move
  # Each move from E1 to E0 is made with probability 1 so that's easy (0 correction on the log scale)

  ### Second type of move
  # Each move from E0 to E1 is made with probability

  # # For the forward move it is:
  # corr_1_E0_to_E1 <- sum(
  #   sapply(seq_along(date_idx_E0_to_E1), function(e) {
  #     compute_p_accept_move_from_E0_to_E1(
  #       i = i,
  #       group_idx = group_idx,
  #       date_idx = date_idx_E0_to_E1[e],
  #       curr_aug_dat = proposed_aug_dat_step1,
  #       proposed_aug_dat = proposed_aug_dat_step2,
  #       theta = theta,
  #       obs_dat = obs_dat,
  #       hyperparameters = hyperparameters,
  #       index_dates = index_dates,
  #       range_dates = range_dates)[2]}))

  # # For the backwards move it is:
  # corr_1_E0_to_E1_rev <-  - sum(
  #   sapply(seq_along(date_idx_E0_to_E1), function(e) {
  #     compute_p_accept_move_from_E0_to_E1(
  #       i = i,
  #       group_idx = group_idx,
  #       date_idx = date_idx_E0_to_E1[e],
  #       curr_aug_dat = proposed_aug_dat_rev2,
  #       proposed_aug_dat = proposed_aug_dat_rev1,
  #       theta = theta,
  #       obs_dat = obs_dat,
  #       hyperparameters = hyperparameters,
  #       index_dates = index_dates,
  #       range_dates = range_dates)[2]}))

  ### Third type of move
  # Moving missing data
  # Calculation of correction factors is made above

  ### Recap
  # probability of moves from E1 to E0 is 1 so corresponding term is zero
  # corr_1_E0_to_E1 # computes (- log_p_proposing_move_E0old_to_E1new)
  # corr_1_E0_to_E1_rev # computes (- log_p_proposing_move_from_E0new_to_E1old)
  # correct_factor_new_delay # computes the log_prob_proposing_new_missing_values
  # correct_factor_new_delay_rev # computes the log_prob_proposing_old_missing_values

  corr <- (corr_1_E0_to_E1 - corr_1_E0_to_E1_rev) + # to obtain log_p_new_to_old - log_p_old_to_new
    (correct_factor_new_delay_rev - correct_factor_new_delay) # to obtain log_p_new_to_old - log_p_old_to_new

  #print(paste(c(ratio_post ,
  #                corr_1_E0_to_E1 , corr_1_E0_to_E1_rev ,
  #              correct_factor_new_delay_rev , correct_factor_new_delay)))

  # if(is.infinite(ratio_post) & is.infinite(corr)) {
  #   browser()
  #   # this seems to yield problems currently when two dates are wrong and one right
  # }

  p_accept <- ratio_post + corr

  if (p_accept > 0) p_accept <- 0
  # print(p_accept)

  # accept/reject step
  tmp <- log(runif(1))
  if (tmp < p_accept) { # accepting with a certain probability
    new_aug_dat <- proposed_aug_dat_step3
    accept <- 1
  } else { # reject
    new_aug_dat <- curr_aug_dat
    accept <- 0
  }

  # return a list of size 2 where
  #		the first value is the new augmented data set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  res <- list(new_aug_dat = new_aug_dat, accept = accept)

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
#' \itemize{
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
#' \itemize{
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
#'  \itemize{
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
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param hyperparameters A list of hyperparameters: see details.
#' @details \code{curr_theta} should be a list containing:
#' \itemize{
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
#' \itemize{
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
#' \itemize{
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
#'  \itemize{
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
                            obs_dat,
                            hyperparameters) {

  tmp <- compute_n_errors(aug_dat, obs_dat)
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
# test_move_zeta_gibbs <- move_zeta_gibbs(aug_dat, curr_theta = theta, obs_dat,
# hyperparameters)
# test_move_zeta_gibbs$new_theta$zeta # new value
# theta$zeta # old value
