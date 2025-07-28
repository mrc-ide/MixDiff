#######################################
### Functions to simulate a dataset ###
#######################################

#' Simulates a dataset with true event dates, and optionally adds missingness
#'  and errors.
#' 
#' @param theta List of parameters for delay distributions and data errors (see
#'  details).
#' @param n_per_group Vector containing the number of individuals to simulate
#'  per group (i.e., community-alive, hospitalised-dead, etc.).
#' @param range_dates Vector of two integers: the range to draw the initial date
#'  from (uniformly).
#' @param index_dates A list defining how each group's date delays are simulated
#'  (see details).
#' @param simul_error Boolean. If TRUE, simulate missing and erroneous data.
#' @param remove_allNA_indiv Boolean. If TRUE, remove rows for individuals with
#'  all dates missing (only used if \code{simul_error = TRUE}).
#'  
#' @return A list with three elements: 
#'  \itemize{
#'    \item{\code{true_dat}: List of matrices containing the simulated true
#'     dates as integers for each group. Each row corresponds to an individual;
#'      columns correspond to event dates (e.g., onset, report) as defined in
#'       \code{index_dates}.}
#'    \item{\code{obs_dat}: Same structure as \code{true_dat} but with simulated
#'     errors and missing dates (NULL if \code{simul_error = FALSE}).}
#'     \item{\code{E}: Error indicator matrices with the same structure as
#'      \code{true_dat} and \code{obs_dat}, where each element indicates: -1 for
#'       missing date, 1 for erroneous data, or 0 for correctly recorded data. E
#'        will be NULL if \code{simul_error = FALSE}.}
#'  }
#'  
#' @details
#' \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}: A list of length \code{n_groups} (the number of groups
#'   to be simulated). Each element should be a scalar or vector giving
#'    the mean delay(s) used to simulate dates in that group.}
#'  \item{\code{CV}: A list of length \code{n_groups}. Each element of \code{CV}
#'   should be a scalar or vector giving the coefficient of variation for the
#'    delay(s) used to simulate dates in that group.}
#'  \item{\code{prop_missing_data} (only required if \code{simul_error = TRUE}):
#'   A scalar in [0,1] giving the probability that a date is missing.}
#'  \item{\code{zeta} (only required if \code{simul_error = TRUE}): A scalar
#'   in [0,1] giving the probability that, if a date is not missing, it is
#'    recorded with error.}
#' }
#' 
#' \code{index_dates} should be a list of length \code{n_groups}. Each element
#'  should be a matrix with 2 rows and columns corresponding to the delays of
#'   interest for that group. For each column:
#'   \itemize{
#'    \item Row 1 gives the index of the origin date.
#'    \item Row 2 gives the index of the destination date.
#'    }
#' The number of columns in \code{index_dates[[k]]} should match the length of
#'  \code{theta$mu[[k]]} and \code{theta$CV[[k]]}.
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and
#'  c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are
#'   respectively the mean and coefficient of variation of two delays: the
#'    first delay being between date 1 and date 2, and the second being between
#'     date 1 and date 3.
#'
#' In the simulation, date 1 will be drawn uniformly within \code{range_dates}. 
#' Then date 2 will be drawn as date 1 + a delay drawn from a discretised gamma
#'  distribution with mean theta$mu[[k]][1] and CV theta$CV[[k]][1]. 
#' Finally, date 3 will be drawn as date 1 + a delay drawn from a discretised
#'  gamma distribution with mean theta$mu[[k]][1] and CV theta$CV[[k]][1].
#'  
#'  @seealso
#'    [simul_obs_dat()] for simulating missing and erroneous data,
#'    [date_to_int()] for converting dates to integers corresponding to the
#'     number of days from a given origin,
#'    [int_to_date()] for converting integers to dates based on a given origin
#'     from which the integer counts the number of days.
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
#' # Set up the parameters for the simulation
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of delays
#' 
#' # Number of individuals to simulate in each group
#' n_per_group <- rep(10, n_groups)
#' 
#' # Range of dates in which to draw first set of dates for each group
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
#'                              as.Date("01/01/2015", "%d/%m/%Y")))
#' 
#' # Delays to use to simulate subsequent dates from the first, in each group
#' index_dates <- list(matrix(c(1, 2), nrow = 2), cbind(c(1, 2), c(1, 3)))
#' 
#' # Perform the simulation
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
simul_true_data <- function(
    theta,
    n_per_group,
    range_dates,
    index_dates,
    simul_error = FALSE,
    remove_allNA_indiv = TRUE
) {

  # Check input
  if (simul_error) {
    required_params <- c("prop_missing_data", "zeta", "mu", "CV")
    missing_params <- setdiff(required_params, names(theta))
    
    if (length(missing_params) > 0) {
      stop(
        sprintf(
          "When `simul_error = TRUE`, `theta` must include the following parameters: %s. Missing: %s.",
          paste(required_params, collapse = ", "),
          paste(missing_params, collapse = ", ")
        )
      )
    }
  }
  
  # Initialise output list for each group
  D <- list()
  
  # Loop through each group
  for (g in seq_along(theta$mu)) {
    
    # Simulate 20% more individuals per group than needed in case of all-NA rows
    extra_rows <- n_per_group[g] * 1.2
    
    # Initialise matrix: cols = number of events per person, rows = individuals
    D[[g]] <- matrix(NA, extra_rows, length(theta$mu[[g]]) + 1)
    
    # Simulate initial event date uniformly from range
    D[[g]][, 1] <- sample(
      seq(range_dates[1], range_dates[2], 1), extra_rows, replace = TRUE
      )

    # Generate all subsequent dates from delays (via gamma distribution)
    for (j in seq_len(ncol(index_dates[[g]]))) {
      mu <- theta$mu[[g]][j]
      CV <- theta$CV[[g]][j]
      
      # Sample from discretised gamma
      delay <- discr_gamma_sample(extra_rows, mu, CV)
      
      # Compute new date = origin date + delay
      D[[g]][, index_dates[[g]][2, j]] <-
        D[[g]][, index_dates[[g]][1, j]] + delay
      
    }
  }
  
  # If simulating errors and/or missing data create an "observed" dataset
  if (simul_error) {
    observed_D <- simul_obs_dat(D,
                                theta,
                                range_dates,
                                remove_allNA_indiv = TRUE,
                                n_per_group)
    
    return(list(true_dat = observed_D$true_dat,
                obs_dat = observed_D$obs_dat,
                E = observed_D$E))
    
  } else {
    
    # Loop through each group
    for (g in seq_along(theta$mu)) {
      
    # Remove excess simulated individuals so that nrow == n_per_group
    # (this is done within simul_obs_dat otherwise)
    if (nrow(D[[g]]) > n_per_group[g]) {
      D[[g]] <- D[[g]][1:n_per_group[g], ]
    }
      }
    
    return(list(true_dat = D, obs_dat = NULL, E = NULL))
    
  }
}

#' Create dataset which introduces missingness and errors.
#' 
#' @param D A list of true data (\code{true_dat}) matrices returned by
#'  \code{\link{simul_true_data}}). 
#' @param theta A list of parameters including:
#' \itemize{
#'  \item{\code{prop_missing_data}: probability of each date being missing}
#'  \item{\code{zeta}: probability that a non-missing date is recorded with error}
#'  }
#' @param range_dates The range of dates (as integers) to draw the erroneous
#'  values from (uniformly).
#' @param remove_allNA_indiv Boolean. If TRUE, individuals with only NA dates
#'  will be removed.
#'
#' @return A list with two elements:
#' \itemize{
#'  \item{\code{obs_dat}: Same structure as \code{D}, but where some dates are
#'   now missing, and some are erroneous.}
#'  \item{\code{E}: Error indicator matrices with the same structure as
#'   \code{D} and \code{obs_dat}, where each element indicates: -1 for a
#'    missing date, 1 for an erroneous date, or 0 for a correctly recorded
#'     date.}
#' }
#'
#' @details
#' \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{prop_missing_data} (only required if \code{simul_error = TRUE}):
#'   A scalar in [0,1] giving the probability that a date is missing.}
#'  \item{\code{zeta} (only required if \code{simul_error = TRUE}): A scalar
#'   in [0,1] giving the probability that, if a date is not missing, it is
#'    recorded with error.}
#' }
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
#' theta$CV <- list(0.5, c(0.5, 0.5)) # CV of these delays
#' theta$prop_missing_data <- 0.25 # prob of date missing in observations
#' theta$zeta <- 0.05 # prob that non-missing date is recorded with error
#' 
#' # Number of individuals to simulate in each group
#' n_per_group <- rep(10, n_groups)
#' 
#' # Range of dates in which to draw the first set of dates for each group
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
#'                              as.Date("01/01/2015", "%d/%m/%Y")))
#' 
#' # Delays to use to simulate subsequent dates from the first, in each group
#' index_dates <- list(matrix(c(1, 2), nrow = 2), cbind(c(1, 2), c(1, 3)))
#' 
#' # Perform the simulation
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' observed_D <- simul_obs_dat(D$true_dat, theta, range_dates,
#'                             remove_allNA_indiv = TRUE)
simul_obs_dat <- function(D,
                          theta,
                          range_dates,
                          remove_allNA_indiv = TRUE,
                          n_per_group) {
  
  # Initialise E to store error indicators (same structure as D)
  E <- D
  
  # Initialise observed dataset as copy of the true data (D)
  obs_dat <- D
  
  # Loop through each group
  for (g in seq_len(length(D))) {
    
    # Loop through each column (each event e.g. onset)
    for (j in seq_len(ncol(D[[g]]))) {
      
      # Randomly assign each value as missing (-1), error (1) or correct (0)
      # based on prop_missing_data and zeta
      E[[g]][,j] <- sample(
        c(-1, 1, 0), nrow(D[[g]]), replace = TRUE,
        prob = c(
          theta$prop_missing_data,
          (1 - theta$prop_missing_data) * theta$zeta,
          (1 - theta$prop_missing_data) * (1 - theta$zeta))
      )
      
      # Set observed value to NA where the date should be missing (E = -1)
      obs_dat[[g]][E[[g]][, j] == -1, j]  <- NA
      
      # Copy true date if correctly observed (E = 0)
      obs_dat[[g]][E[[g]][, j] == 0, j]  <- D[[g]][E[[g]][, j] == 0, j]
      
      # Replace with random date if recorded with error (E = 1)
      # This will need updating if error model changes
      
      # obs_dat[[g]][E[[g]][, j] == 1, j]  <- sample(
      #   seq(range_dates[1], range_dates[2], 1),
      #   sum(E[[g]][, j] == 1),
      #   replace = TRUE
      #   )
      
      # Get the row indices where errors are to be added
      err_idx <- which(E[[g]][, j] == 1)
      true_vals <- D[[g]][err_idx, j]
      range_pool <- seq(range_dates[1], range_dates[2], 1)
      
      err_vals <- integer(length(err_idx))
      
      # Check that the sampled date doesn't match true date
      for (i in seq_along(err_idx)) {
        repeat {
          candidate <- sample(range_pool, 1)
          if (candidate != true_vals[i]) {
            err_vals[i] <- candidate
            break
          }
        }
      }
      
      obs_dat[[g]][err_idx, j] <- err_vals
      
    }
    
    # Remove individuals with all missing dates
    if (remove_allNA_indiv) {
      exclude <- which(rowSums(is.na(obs_dat[[g]])) == ncol(obs_dat[[g]]))
      
      if (length(exclude) > 0) {
        obs_dat[[g]] <- obs_dat[[g]][-exclude, ]
        E[[g]] <- E[[g]][-exclude, ]
        D[[g]] <- D[[g]][-exclude, ]
      }
    }
    
    # Remove excess simulated individuals so that nrow == n_per_group
    if (nrow(obs_dat[[g]]) > n_per_group[g]) {
      obs_dat[[g]] <- obs_dat[[g]][1:n_per_group[g], ]
      E[[g]] <- E[[g]][1:n_per_group[g], ]
      D[[g]] <- D[[g]][1:n_per_group[g], ]
    }
  }
  
  return(list(true_dat = D, obs_dat = obs_dat, E = E))
}
