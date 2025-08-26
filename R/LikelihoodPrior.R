#------------------------------------------------------------------------------
# likelihood function
#------------------------------------------------------------------------------

#' Compute log-likelihood of observed dates for a subset of individuals and dates
#'
#' @description Calculates the likelihood of the observed dates (\code{obs_dat})
#'  given the augmented dates (\code{aug_dat}) for a specific group, set of
#'   individuals, and dates.
#'
#' @details Assumes that:
#'  - if no error (E = 0), true date must equal observed date.
#'  - If error or missing (E = 1 or -1), observed date is treated as uniformly
#'   likely over \code{range_dates}.
LL_observation_term_by_group_delay_and_indiv <- function(aug_dat,
                                                         obs_dat,
                                                         group_idx,
                                                         date_idx,
                                                         indiv_idx,
                                                         range_dates = NULL) {

  if (is.null(range_dates)) range_dates <- find_range(obs_dat)

  LL <- matrix(NA, length(indiv_idx), length(date_idx))

  # Identify where there are no recording errors (E = 0)
  indicator_no_error <- aug_dat$E[[group_idx]][indiv_idx, date_idx] == 0
  no_error <- which(indicator_no_error, arr.ind = TRUE)

  # Where there is not an error the date in aug_dat should match the observed
  # date. If match: log(1) = 0, if mismatch: log(0) = -Inf
  LL[no_error] <- log(aug_dat$D[[group_idx]][indiv_idx, date_idx][no_error] ==
                        obs_dat[[group_idx]][indiv_idx, date_idx][no_error])

  # Handle errors (E = 1) or missingness (E = -1)
  error_or_missing <- which(!indicator_no_error)

  # Likelihood for erroneous/missing data assumed uniform over date range.
  # K = rel prob of observing a given error, conditional on presence of error.
  # For now, K is given as 1/n, where n is the number of dates in the
  # range_dates. Could use something different if we define the space of
  # possible errors differently. Think about impact of this choice.
  K <- (1 / as.numeric(diff(range_dates)))
  dates_to_check <- aug_dat$D[[group_idx]][indiv_idx, date_idx][error_or_missing]
  in_range <- dates_to_check >= range_dates[1] & dates_to_check <= range_dates[2]

  LL[error_or_missing] <- log(K * in_range)

  # Prevent -Inf in posterior by replacing with arbitrary large negative value
  LL[is.infinite(LL)] <- -1e5

  return(LL)
}

#' Compute total observation likelihood over all individuals, dates, and groups
LL_observation_term <- function(aug_dat, theta, obs_dat, range_dates = NULL) {

  if (is.null(range_dates)) range_dates <- find_range(obs_dat)

  LL <- sum(unlist(lapply(
    seq_len(length(obs_dat)),
    function(g) {
      sum(LL_observation_term_by_group_delay_and_indiv(
        aug_dat, theta, obs_dat, g,
        seq_len(ncol(aug_dat$D[[g]])),
        seq_len(nrow(obs_dat[[g]])),
        range_dates
      ))
      }
    )))
  return(LL)
}

#' Compute log-likelihood for recording errors
#'
#' @description
#' For a given group, date index and set of individuals, compute the
#' log-likelihood of the observed error indicators (E). This reflects the
#'  probability of observing a date entry recorded without error, given
#'   the error rate parameter, zeta.
#'
#' @param aug_data List containing augmented data, including:
#' - D: date matrices (as integers) for each group
#' - E: error indicator matrices for each group
#' @param theta List of model parameters, including zeta.
#' @param group_idx Index of the group to consider
#' @param date_idx Indices of the date column(s) to consider
#' @param indiv_idx Indices of the individual(s) to consider
#'
#' @return Matrix of log-likelihood contributions (rows: individuals,
#'  cols: dates)
#'
LL_error_term_by_group_delay_and_indiv <- function(aug_dat,
                                                   theta,
                                                   group_idx,
                                                   date_idx,
                                                   indiv_idx) {

  # Initialise log-likelihood matrix with 0s
  res <- matrix(0, length(indiv_idx), length(date_idx))

  # Identify which entries are not missing (E != -1)
  non_missing <- which(aug_dat$E[[group_idx]][indiv_idx, date_idx] != -1,
                       arr.ind = TRUE)

  # Extract error indicators at those positions (0 = correct, 1 = error)
  tmp <- aug_dat$E[[group_idx]][indiv_idx, date_idx][non_missing]

  # Compute log-likelihood
  # Prob date recorded incorrectly = zeta
  # Prob date recorded correctly = 1 - zeta
  res[non_missing] <- log(theta$zeta) * tmp + log(1 - theta$zeta) * (1 - tmp)

  return(res)
}

#' Compute number of errors and number of observed entries
#' @description
#' Computes the total number of errors and non-missing dates across all
#'  individuals, groups and time points in the augmented dataset.
#'
#' @param aug_dat List containing \code{E}, which is a list of matrices. Each
#'  matrix corresponds to a group and contains the error indicators for each
#'   date per individual (1, 0, -1).
#'
#' @returns Numeric vector of length 2 containing the number of errors and the
#'  number of recorded dates in \code{aug_dat}.
#'
#' @export
#'
#' @examples
#' E_list <- list(matrix(c(1, 0, -1, 1, 0, 1), nrow = 2))
#' aug_dat <- list(E = E_list)
#' compute_n_errors(aug_dat, NULL)
compute_n_errors <- function(aug_dat) {
  number_of_errors <- sum(unlist(aug_dat$E) == 1)
  number_of_recorded_dates <- sum(unlist(aug_dat$E) != -1)
  return(c(number_of_errors, number_of_recorded_dates))
}

#' Compute total log-likelihood for observed errors
LL_error_term <- function(aug_dat, theta) {
  tmp <- compute_n_errors(aug_dat)
  number_of_errors <- tmp[1]
  number_of_recorded_dates <- tmp[2]

  result <- log(theta$zeta) * number_of_errors +
    log(1 - theta$zeta) * (number_of_recorded_dates - number_of_errors)

  return(result)
}

#' Compute delay likelihood for individual delay observations
LL_delays_term_by_group_delay_and_indiv <- function(aug_dat,
                                                    theta,
                                                    obs_dat,
                                                    group_idx,
                                                    delay_idx,
                                                    indiv_idx,
                                                    index_dates,
                                                    Delta = NULL) {

  if (is.null(Delta)) {
    Delta <- compute_delta_group_delay_and_indiv(
      aug_dat$D, group_idx, indiv_idx, delay_idx, index_dates
    )
  }

  # Delay likelihood via discretised gamma
  LL <- DiscrGamma(Delta,
                   mu = theta$mu[[group_idx]][delay_idx],
                   cv = theta$CV[[group_idx]][delay_idx],
                   log = TRUE)
  return(LL)
}

#' Compute total delay likelihood across all individuals and groups
LL_delays_term <- function(aug_dat, theta, obs_dat, index_dates, Delta = NULL) {

  if (is.null(Delta)) {
    Delta <- compute_delta(aug_dat$D, index_dates)
  }

  LL <- sum(sapply(
    seq_along(obs_dat), function(g) { # ANNE: this loops over the groups
      sum(sapply(
        ##ERROR##seq(2, ncol(aug_dat$D[[g]])), # ANNE: this loops over the dates for this group- question: why does it start at 2??
        seq(1, ncol(index_dates[[g]])),
        function(j) {
          sum(
            # ANNE: this is the structure of arguments for LL_delays_term_by_group_delay_and_indiv
            # (aug_dat,theta,obs_dat,group_idx,delay_idx,
            # indiv_idx,index_dates,Delta = NULL)
            LL_delays_term_by_group_delay_and_indiv(
              ##ERROR##aug_dat, theta, obs_dat, g, j - 1,
              aug_dat, theta, obs_dat, g, j,
              seq_len(nrow(obs_dat[[g]])),
              index_dates,
              ##ERROR##Delta[[g]][, j - 1]
              Delta[[g]][, j]
            )
          )
        }
      ))
    }
  ))
  return(LL)
}

#' Compute full likelihood (observation + delay + error)
LL_total <- function(aug_dat, theta, obs_dat, index_dates, range_dates = NULL) {

  res <- LL_observation_term(aug_dat, theta, obs_dat, range_dates) +
    LL_error_term(aug_dat, theta) +
    LL_delays_term(aug_dat, theta, obs_dat, index_dates)

  return(res)
}


#------------------------------------------------------------------------------
# Priors
#------------------------------------------------------------------------------

#' Log prior for error probability parameter zeta (Beta prior)
lprior_prob_error <- function(theta, hyperparameters) {
  dbeta(theta$zeta,
        hyperparameters$shape1_prob_error,
        hyperparameters$shape2_prob_error,
        log = TRUE)

}

#' Log prior for mu or CV (Exponential prior)
lprior_params_delay <- function(what = c("mu", "CV"), theta, hyperparameters) {
  what <- match.arg(what)
  sum(dexp(unlist(theta[[what]]),
           rate = 1 / hyperparameters$mean_mean_delay,
           log = TRUE))
}

#' Total log prior (zeta, mu, CV)
lprior_total <- function(theta, hyperparameters) {
  lprior_prob_error(theta, hyperparameters) +
    lprior_params_delay("mu", theta, hyperparameters) +
    lprior_params_delay("CV", theta, hyperparameters)
}


#------------------------------------------------------------------------------
# Posterior
#------------------------------------------------------------------------------

#' Compute the log joint posterior distribution (likelihood + prior) of
#'  augmented data and parameters given observed data
#'
#' @param aug_dat A list of augmented data with dates `D` and error indicators
#'  `E` for each group.
#' @param theta Parameter list containing: mu, CV, zeta.
#' @param obs_dat A list of observed data in the same format as \code{aug_dat}.
#' @param hyperparameters List of priors for mu, CV, zeta.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}.
#'  If NULL, will be computed automatically.
#'
#' @details
#' \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}: A list of length \code{n_groups}. Each element of \code{mu}
#'   should be a scalar or vector giving the mean delay(s) to use for the
#'    simulation of dates in that group.}
#'  \item{\code{CV}: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar or vector giving the coefficient of variation
#'    of the delay(s) to use to simulate dates in that group.}
#'  \item{\code{zeta}: A scalar in [0;1] giving the probability that, if a
#'   date is not missing, it is recorded with error.}
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
#' @return A scalar giving the value of the log posterior.
#'
#' @export
#'
#' @examples
#' # Number of groups of individuals to simulate
#' n_groups <- 2
#'
#' # Number of dates to simulate for each group
#' n_dates <- c(2, 3)
#'
#' # Setting up the parameters for the simulation
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of delays
#' theta$prop_missing_data <- 0.25 # probability data missing in observations
#' theta$zeta <- 0.05 # probability that non-missing date is recorded with error
#'
#' # Number of individuals to simulate in each group
#' n_per_group <- rep(10, n_groups)
#'
#' # Range of dates in which to draw the first set of dates for each group
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
#'                              as.Date("01/01/2015", "%d/%m/%Y")))
#'
#' # Delays used to simulate subsequent dates from the first, in each group
#' index_dates <- list(matrix(c(1, 2), nrow = 2), cbind(c(1, 2), c(1, 3)))
#'
#' # Simulate data
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' observed_D <- simul_obs_dat(D$true_dat, theta, range_dates,
#'                             remove_allNA_indiv = TRUE)
#' obs_dat <- observed_D$obs_dat
#' true_aug_dat <- list(D = D$true_dat, E = observed_D$E)
#'
#' # Define hyperparameters
#' hyperparameters <- list(shape1_prob_error = 3, shape2_prob_error = 12,
#'                      mean_mean_delay = 100, mean_CV_delay = 100)
#'
#' # Compute log posterior distribution for that data
#' lposterior_total(true_aug_dat, theta, obs_dat, hyperparameters, index_dates,
#'                  range_dates = NULL)
#'
#' # Now use initialised augmented data and check that posterior value for this
#' # is lower than for true data:
#' MCMC_settings <- list(init_options = list(mindelay = 0, maxdelay = 100))
#' aug_dat <- initialise_aug_data(observed_D$obs_dat, index_dates,
#'                                MCMC_settings)
#' lposterior_total(aug_dat, theta, obs_dat, hyperparameters, index_dates,
#'                  range_dates = NULL)
lposterior_total <- function(aug_dat, theta, obs_dat, hyperparameters,
                             index_dates, range_dates = NULL) {

  LL_total(aug_dat, theta, obs_dat, index_dates, range_dates) +
    lprior_total(theta, hyperparameters)
}
