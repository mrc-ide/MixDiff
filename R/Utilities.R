###----------------------------------------------------------------------------
### discrete Gamma Probability Mass Function - see Cori et al. AJE 2013
###----------------------------------------------------------------------------

#' Computes the probability mass of the discrete gamma distribution at integer
#'  k.
#'
#' @param k Integer (or vector of integers).
#' @param mu Mean of the gamma distribution.
#' @param cv Coefficient of variation.
#' @param sigma Optional standard deviation. Defaults to mu * cv.
#' @param log Logical; if TRUE, returns log-probabilities.
#'
#' @return Probability or log-probability mass at k.
#' @export
DiscrGamma <- function(k, mu, cv = NULL, sigma = mu * cv, log = TRUE) {
  
  if (!is.null(cv)) {
    if (cv < 0) stop("cv must be >= 0.")
  }
  if (sigma < 0) stop("sigma must be >= 0.")

  shape <- (mu / sigma)^2
  rate <- mu / (sigma^2)

  res <- (k + 1) * pgamma(k + 1, shape, rate) +
    (k - 1) * pgamma(k - 1, shape, rate) -
    2 * k * pgamma(k, shape, rate)

  res <- res + (shape / rate) * (
    2 * pgamma(k, shape + 1, rate) -
      pgamma(k - 1, shape + 1, rate) -
      pgamma(k + 1, shape + 1, rate)
  )

  res <- pmax(0, res)
  if (log) return(log(res)) else return(res)
}

# Copied across from likelihood prior in case helpful
# DiscrGamma_approx <- function (k, mu, CV = NULL, sigma=mu*CV, log=TRUE) 
# {
#   if (!is.null(CV)) {
#     if(CV < 0)
#       stop("CV must be >=0.")
#   }
#   if (sigma < 0) {
#     stop("sigma must be >=0.")
#   }
#   shape <- (mu/sigma)^2
#   rate <- mu/(sigma^2)
#   
#   res <- distcrete::distcrete("gamma", 1, shape, rate, w = 0.5)$d(k)
#   return(if(log) log(res) else res)
# }


#' Sample from the Discrete Gamma Distribution
#'
#' @param n Number of samples.
#' @param mu Mean of the gamma.
#' @param cv Coefficient of variation.
#'
#' @return Integer vector of samples.
#' @export
discr_gamma_sample <- function(n, mu, cv) {
  
  # Convert mean and cv to gamma parameters
  shape <- (mu / (mu * cv))^2
  rate <- mu / (mu * cv)^2
  
  # 99.9th percentile of the continuous gamma distribution (max day)
  k_max <- ceiling(qgamma(0.999, shape = shape, rate = rate))
  ks <- 0:k_max
  
  # Compute the discrete probability mass for each k
  probs <- pmax(0, DiscrGamma(ks, mu, cv, log = FALSE))
  
  # Normalise to sum to 1
  probs <- probs / sum(probs)

  # Sample from the discrete distribution
  sample(ks, size = n, replace = TRUE, prob = probs)
}


###----------------------------------------------------------------------------
### functions to handle dates
###----------------------------------------------------------------------------

#' Convert date to integer corresponding to number of days from a given origin
#'
#' @param date \code{Date} object to be converted to integer.
#' @param origin The date used as origin for counting days.
#' @return An integer corresponding to the number of days since \code{origin}.
#' @export
#' @examples
#' date_to_int(as.Date("1970-01-01"))
date_to_int <- function(date, origin = "1970-01-01") {
  return(as.integer(date - as.Date(origin)))
}

#' Convert integer to date, based on a given origin from which the integer
#'  counts the number of days
#'
#' @param int integer to be converted to \code{Date}.
#' @param origin The date used as origin for counting days.
#' @return A \code{Date} object corresponding to \code{int} days after
#'  \code{origin}.
#' @export
#' @examples
#' int_to_date(365)
int_to_date <- function(int, origin = "1970-01-01") {
  return(int + as.Date(origin))
}

###----------------------------------------------------------------------------
### functions to find parameters of distributions from mean/var or mean/sd
###----------------------------------------------------------------------------

#' Finds the two shape parameters of the beta distribution, given a mean and
#'  variance
#'
#' @param mean Mean of the beta distribution
#' @param var Variance of the beta distribution
#' @return A vector containing the two shape parameters of the beta
#'  distribution.
#' @export
#' @examples
#' param_beta <- find_params_beta(0.1, 0.05)
#' sample_beta <- rbeta(1000, param_beta[1], param_beta[2])
#' mean(sample_beta) # compare to 0.1
#' var(sample_beta) # compare to 0.05
find_params_beta <- function(mean, var) {
  # for a beta distribution:
  # mean = shape1 / (shape1 + shape2)
  # var = shape1 * shape2 / ((shape1 + shape2)^2 * (shape1 + shape2 + 1))
  # after solving we find that
  shape1 <- mean^2 * (1 - mean) / var - mean
  shape2 <- shape1 * (1 / mean - 1)
  return(c(shape1, shape2))
}

#' Finds the shape and scale parameters of the gamma distribution, given a mean
#'  and standard deviation
#'
#' @param mean Mean of the gamma distribution
#' @param sigma Standard deviation of the gamma distribution
#' @param CV An alternative way to specify the stndard deviation, through the
#'  coefficient of variation, that is \code{sigma}/\code{mean}
#' @return A vector containing the shape and scale parameters of the gamma
#'  distribution.
#' @export
#' @examples
#' param_gamma <- find_params_gamma(mean = 0.1, sigma = 0.05)
#' sample_gamma <- rgamma(1000, shape = param_gamma[1], scale = param_gamma[2])
#' mean(sample_gamma) # compare to 0.1
#' sd(sample_gamma) # compare to 0.05
find_params_gamma <- function(mean, sigma = mean * CV, CV) {
  shape <- (mean / sigma)^2
  scale <- sigma^2 / mean
  return(c(shape, scale))
}

###----------------------------------------------------------------------------
### compute_delta functions to compute relevant delays based on index, which
### tells you which dates should be used for delayl calculation
###----------------------------------------------------------------------------

#' Compute relevant delays for one individual, based on data
#' 
#' @param D A list of data, in the format of the first element (called
#'  \code{true_dat}) in the list returned by \code{\link{simul_true_data}}. 
#' @param group_idx a scalar or vector giving the index of the group(s) for
#'  which to calculate the delays
#' @param indiv_idx a scalar or vector giving the index of the individuals in
#'  group \code{group_idx} for which to calculate the delays
#' @param delay_idx a scalar or vector giving the index of the delay(s) to be
#'  calculated
#' @param index_dates A list containing indications on which delays to consider
#'  in the simulation. 
#' @details \code{index_dates} should be a list of length
#'  \code{n_groups=length(D)}. Each element of \code{index_dates} should be a
#'   matrix with 2 rows and a number of columns corresponding to the delays of
#'    interest for that group. 
#' For each column (i.e. each delay), the first row gives the index of the
#'  origin date, and the second row gives the index of the destination date. 
#' @return A list of same length as \code{D}. Each element in the list is a
#'  matrix with same number of rows as in \code{D}, but potentially a different
#'   number of columns, corresponding to the relevant delays calculated for that
#'    group. 
#' @export
#' @examples
#' # Number of groups of individuals to simulate
#' n_groups <- 3
#' # Number of dates to simulate for each group
#' n_dates <- c(2, 3, 4)
#' # Setting up the parameters for the simulation
#' theta <- list()
#' theta$mu <- list(5, c(10, 15), c(10, 15, 20)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5)) # cv of these delays
#' # Number of individuals to simulate in each group
#' n_per_group <- rep(10, n_groups)
#' # Range of dates in which to draw the first set of dates for each group
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
#'                              as.Date("01/01/2015", "%d/%m/%Y")))
#' # Which delays to use to simulate subsequent dates from the first, in each
#'  # group?
#' index_dates <- list(matrix(c(1, 2), nrow = 2),
#'                     cbind(c(1, 2), c(1, 3)),
#'                     cbind(c(1, 2), c(2, 3),c(1, 4)))
#' # Perform the simulation
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' # Compute the first delay for first individual in first group
#' compute_delta_group_delay_and_indiv(D$true_dat, group_idx = 3, indiv_idx = 1,
#'  delay_idx = 3, index_dates)
#'
compute_delta_group_delay_and_indiv <- function(D, group_idx, indiv_idx,
                                                delay_idx, index_dates) {
  
  Delta <- D[[group_idx]][indiv_idx, index_dates[[group_idx]][, delay_idx][2]] -
    D[[group_idx]][indiv_idx, index_dates[[group_idx]][, delay_idx][1]]
  
  return(Delta)
}

#' Compute relevant delays based on data
#' 
#' @param D A list of data, in the format of the first element (called
#'  \code{true_dat}) in the list returned by \code{\link{simul_true_data}}. 
#' @param index_dates A list containing indications on which delays to consider
#'  in the simulation. 
#' @details \code{index_dates} should be a list of length
#'  \code{n_groups = length(D)}. Each element of \code{index_dates} should be a
#'   matrix with 2 rows and a number of columns corresponding to the delays of
#'    interest for that group. 
#' For each column (i.e. each delay), the first row gives the index of the
#'  origin date, and the second row gives the index of the destination date. 
#' @return A list of same length as \code{D}. Each element in the list is a
#'  matrix with same number of rows as in \code{D}, but potentially a different
#'   number of columns, corresponding to the relevant delays calculated for that
#'    group. 
#' @export
#' @examples
#' ### Number of groups of individuals to simulate
#' n_groups <- 3
#' ### Number of dates to simulate for each group
#' n_dates <- c(2, 3, 4)
#' ### Setting up the parameters for the simulation
#' theta <- list()
#' theta$mu <- list(5, c(10, 15), c(10, 15, 20)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5)) # cv of these delays
#' ### Number of individuals to simulate in each group
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
#'  as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each
#'  # group?
#' index_dates <- list(matrix(c(1, 2), nrow = 2),
#'                     cbind(c(1, 2), c(1, 3)),
#'                     cbind(c(1, 2), c(2, 3),c(1, 4)))
#' ### Perform the simulation
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' ### Compute the corresponding delays
#' Delays <- compute_delta(D$true_dat, index_dates)
#'
compute_delta <- function(D, index_dates) {
  
  Delta <- lapply(seq_len(length(D)), function(g) {
    m <- matrix(NA, nrow(D[[g]]), ncol(D[[g]]) - 1)
    for (j in seq_len(ncol(m))) {
      m[, j] <- D[[g]][, index_dates[[g]][, j][2]] -
        D[[g]][, index_dates[[g]][, j][1]]
    }
    return(m)
  })
  return(Delta)
}

###----------------------------------------------------------------------------
### find range of dates
###----------------------------------------------------------------------------

#' Find range of dates from a dataset
#' 
#' @param obs_dat A list of data, in the format of the first element (called
#'  \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @return A vector of two integers coresponding to the range in \code{obs_dat}. 
#' @export
#' @examples
#' # Number of groups of individuals to simulate
#' n_groups <- 2
#' # Number of dates to simulate for each group
#' n_dates <- c(2, 3)
#' # Setting up the parameters for the simulation
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of these delays
#' # Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' # Range of dates in which to draw the first set of dates for each group
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"),
#'  as.Date("01/01/2015", "%d/%m/%Y")))
#' # Which delays to use to simulate subsequent dates from the first, in each
#'  group?
#' index_dates <- list(matrix(c(1, 2), nrow = 2), cbind(c(1, 2), c(1, 3)))
#' # Perform the simulation
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' # Find the range
#' find_range(D$true_dat)
#' # Compare with range specified in the first place
#' range_dates
find_range <- function(obs_dat) {
  
  min_date <- Inf
  max_date <- -Inf
  
  for (g in seq_along(obs_dat)) {
    for (j in seq_len(ncol(obs_dat[[g]]))) {
      
      current_col <- obs_dat[[g]][, j]
      
      if (all(is.na(current_col))) {
        stop(paste0("All dates in column ", j, ", group ", g, " are NA. ",
                    "Check that each date column in obs_dat has >=1 date."))
      }
      
      min_date_tmp <- min(current_col, na.rm = TRUE)
      min_date <- min(c(min_date, min_date_tmp), na.rm = TRUE)
      
      max_date_tmp <- max(current_col, na.rm = TRUE)
      max_date <- max(c(max_date, max_date_tmp), na.rm = TRUE)
    }
  }
  return(c(min_date, max_date))
}


###----------------------------------------------------------------------------
### compute rules on the order of dates based on index_dates object
###----------------------------------------------------------------------------

#' Compute rules on the order of dates based on index_dates object; see details
#' 
#' @param index_dates A list containing indications on which delays to consider
#'  in the simulation, see details.
#' @details \code{index_dates} should be a list; each elements corresponding to
#'  a group of individuals of interest. Each element of \code{index_dates}
#'   should be a matrix with 2 rows and a number of columns corresponding to
#'    the delays of interest for that group. For each column (i.e. each delay),
#'     the first row gives the index of the origin date, and the second row
#'      gives the index of the destination date. 
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and
#'  c(1, 3), this indicates that for group \code{k} we are interested in two
#'   delays: the first delay being between date 1 and date 2, and the second
#'    being between date 1 and date 3. 
#' 
#' This function is used to find appropriate starting points for the MCMC,
#'  i.e. when choosing initial values for the missing dates, this allows making
#'   sure the chosen value is consistent with the ordering of dates in that
#'    group, and hence will not generate a null likelihood. 
#' 
#' @return A list of same length as index_dates, containing indications on
#'  ordering of dates for each group. More specifically, each element of
#'   \code{index_dates_order} is a matrix with 2 rows and a number of columns
#'    corresponding to the delays with order rules for that group.
#' For each column (i.e. each delay), the first row gives the index of the
#'  origin date, and the second row gives the index of the destination date.
#' Each column specifies a rule saying that the origin date must be before the
#'  destination date.  
#' This is computed from \code{index_dates} assuming all delays specified in
#'  \code{index_dates} have to be positive, and using transitivity rules to
#'   derive potential additional rules of positivity. 
#' In the example below, in group 3, \code{index_dates} indicates that the
#'  delay between dates 1 and 2, and the delay between dates 2 and 3, are
#'   positive. Hence by transitivity, the delay between date 1 and 3 has to be
#'    positive as well, as seen in the output of the function in that example. 
#' @export
#' @examples
#' index_dates <- list(matrix(c(1, 2), nrow = 2), 
#'                        cbind(c(1, 2), c(1, 3)), 
#'                        cbind(c(1, 2), c(2, 3), c(1, 4)), 
#'                        cbind(c(1, 2), c(2, 3), c(1, 4)))
#' index_dates_order <- compute_index_dates_order(index_dates)
compute_index_dates_order <- function(index_dates) {
  
  index_dates_order <- index_dates
  
  # Remove duplicates
  for (e in seq_along(index_dates_order)) {
    index_dates_order[[e]] <- unique(t(index_dates_order[[e]]), MARGIN = 1)
    index_dates_order[[e]] <- t(index_dates_order[[e]])
  }
  
  has_changed <- TRUE
  while (has_changed) {
    has_changed <- FALSE
    for (e in seq_along(index_dates_order)) {
      current_rules <- as.matrix(index_dates_order[[e]])
      
      all_pairs <- t(current_rules)
      new_pairs_list <- list()
      
      # Find new transitive relationships e.g. if A -> B and B -> C then A -> C
      for (i in 1:nrow(all_pairs)) {
        for (j in 1:nrow(all_pairs)) {
          if (all_pairs[i, 2] == all_pairs[j, 1]) {
            new_rule <- c(all_pairs[i, 1], all_pairs[j, 2])
            # Add the new rule to a temp list
            new_pairs_list[[length(new_pairs_list) + 1]] <- new_rule
          }
        }
      }
      
      # Add new rules if any were found
      if (length(new_pairs_list) > 0) {
        # Convert the list of new pairs into a matrix
        new_matrix <- do.call(cbind, new_pairs_list)
        # Combine the original and new rules and then find unique columns
        combined_rules <- cbind(current_rules, new_matrix)
        unique_rules <- unique(t(combined_rules), MARGIN = 1)
        
        # Check if new rules were added
        if (nrow(unique_rules) > ncol(current_rules)) {
          index_dates_order[[e]] <- t(unique_rules)
          has_changed <- TRUE
        }
      }
    }
  }
  
  return(index_dates_order)
}


#' Check MCMC settings
#' 
#' @export
check_MCMC_settings <- function(MCMC_settings, index_dates) {
  if (MCMC_settings$chain_properties$n_iter < MCMC_settings$chain_properties$burnin)
    stop("Burnin must be <= n_iter")
  if (length(MCMC_settings$moves_options$sdlog_mu) != length(index_dates))
    stop("sdlog_mu does not have the correct length")
  if (length(MCMC_settings$moves_options$sdlog_CV) != length(index_dates))
    stop("sdlog_CV does not have the correct length")
  if (!all(lengths(MCMC_settings$moves_options$sdlog_mu) == lengths(index_dates) / 2))
    stop("sdlog_mu does not have the correct structure")
  if (!all(lengths(MCMC_settings$moves_options$sdlog_CV) == lengths(index_dates) / 2))
    stop("sdlog_CV does not have the correct structure")
}
