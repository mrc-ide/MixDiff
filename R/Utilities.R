###############################################
### functions to handle dates ###
###############################################

#' Convert date to integer corresponding to number of days from a given origin
#'
#' @param date \code{Date} object to be converted to integer.
#' @param origin The date used as origin for counting days.
#' @return An integer corresponding to the number of days since \code{origin}.
#' @export
#' @examples
#' date_to_int(as.Date("1970-01-01"))
date_to_int <- function(date, origin = "1970-01-01")
{
  return(as.integer(date - as.Date(origin)))
}

#' Convert integer to date, based on a given origin from which the integer counts the number of days
#' 
#' @param int integer to be converted to \code{Date}.
#' @param origin The date used as origin for counting days.
#' @return A \code{Date} object corresponding to \code{int} days after \code{origin}.
#' @export
#' @examples
#' int_to_date(365)
int_to_date <- function(int, origin = "1970-01-01")
{
  return(int + as.Date(origin))
}

###############################################
### functions to find parameters of distributions from mean/var or mean/sd ###
###############################################

#' Finds the two shape parameters of the beta distribution, given a mean and variance
#' 
#' @param mean Mean of the beta distribution
#' @param var Variance of the beta distribution
#' @return A vector containing the two shape parameters of the beta distribution.
#' @export
#' @examples
#' param_beta <- find_params_beta(0.1, 0.05)
#' sample_beta <- rbeta(1000, param_beta[1], param_beta[2])
#' mean(sample_beta) # compare to 0.1
#' var(sample_beta) # compare to 0.05
find_params_beta <- function(mean, var) # function to determine parameters of the beta distribution corresponding to a given mean and variance
{
  # for a beta distribution:
  # mean = shape1/(shape1+shape2) 
  # var = shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1)) 
  # after solving we find that
  shape1 <- mean^2*(1-mean)/var - mean
  shape2 <- shape1*(1/mean-1)
  return(c(shape1, shape2))
}

#' Finds the shape and scale parameters of the gamma distribution, given a mean and standard deviation
#' 
#' @param mean Mean of the gamma distribution
#' @param sigma Standard deviation of the gamma distribution
#' @param CV An alternative way to specify the stndard deviation, through the coefficient of variation, that is \code{sigma}/\code{mean}
#' @return A vector containing the shape and scale parameters of the gamma distribution.
#' @export
#' @examples
#' param_gamma <- find_params_gamma(mean=0.1, sigma=0.05)
#' sample_gamma <- rgamma(1000, shape=param_gamma[1], scale=param_gamma[2])
#' mean(sample_gamma) # compare to 0.1
#' sd(sample_gamma) # compare to 0.05
find_params_gamma <- function(mean, sigma=mean*CV, CV) # function to determine parameters of the gamma distribution corresponding to a given mean and std
{
  shape <- (mean/sigma)^2
  scale <- sigma^2/(mean)
  return(c(shape, scale))
}

###############################################
### compute_delta functions to compute relevant delays based on index, which tells you which dates should be used for delayl calculation ###
###############################################

#' Compute relevant delays for one individual, based on data
#' 
#' @param D A list of data, in the format of the first element (called \code{true_dat}) in the list returned by \code{\link{simul_true_data}}. 
#' @param group_idx a scalar or vector giving the index of the group(s) for which to calculate the delays
#' @param indiv_idx a scalar or vector giving the index of the individuals in group \code{group_idx} for which to calculate the delays
#' @param delay_idx a scalar or vector giving the index of the delay(s) to be calculated
#' @param index_dates A list containing indications on which delays to consider in the simulation. 
#' @details \code{index_dates} should be a list of length \code{n_groups=length(D)}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. 
#' For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' @return A list of same length as \code{D}. Each element in the list is a matrix with same number of rows as in \code{D}, but potentially a different number of columns, corresponding to the relevant delays calculated for that group. 
#' @export
#' @examples
#' ### Number of groups of individuals to simulate ###
#' n_groups <- 2
#' ### Number of dates to simulate for each group ###
#' n_dates <- c(2, 3)
#' ### Setting up the parameters for the simulation ###
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of these delays
#' ### Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each group? ###
#' index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)))
#' ### Perform the simulation ###
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' ### Compute the first delay for first individual in first group ###
#' compute_delta_group_delay_and_indiv(D$true_dat, group_idx=1, indiv_idx=1, delay_idx=1, index_dates)
compute_delta_group_delay_and_indiv<-function(D, group_idx, indiv_idx, delay_idx, index_dates)
{
  Delta <- D[[group_idx]][indiv_idx,index_dates[[group_idx]][,delay_idx][2]] - D[[group_idx]][indiv_idx,index_dates[[group_idx]][,delay_idx][1]]
  return(Delta)
}

#' Compute relevant delays based on data
#' 
#' @param D A list of data, in the format of the first element (called \code{true_dat}) in the list returned by \code{\link{simul_true_data}}. 
#' @param index_dates A list containing indications on which delays to consider in the simulation. 
#' @details \code{index_dates} should be a list of length \code{n_groups=length(D)}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. 
#' For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' @return A list of same length as \code{D}. Each element in the list is a matrix with same number of rows as in \code{D}, but potentially a different number of columns, corresponding to the relevant delays calculated for that group. 
#' @export
#' @examples
#' ### Number of groups of individuals to simulate ###
#' n_groups <- 2
#' ### Number of dates to simulate for each group ###
#' n_dates <- c(2, 3)
#' ### Setting up the parameters for the simulation ###
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of these delays
#' ### Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each group? ###
#' index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)))
#' ### Perform the simulation ###
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' ### Compute the corresponding delays ###
#' Delays <- compute_delta(D$true_dat, index_dates)
compute_delta <- function(D, index_dates)
{
  Delta <- lapply(1:length(D), function(g){
    m <- matrix(NA, nrow(D[[g]]), ncol(D[[g]])-1)
    for(j in 1: ncol(m))
    {
      m[,j] <- D[[g]][,index_dates[[g]][,j][2]] - D[[g]][,index_dates[[g]][,j][1]]
    }
    return(m)
  })
  return(Delta)
}

###############################################
### find range of dates
###############################################

#' Find range of dates from a dataset
#' 
#' @param obs_dat A list of data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @return A vector of two integers coresponding to the range in \code{obs_dat}. 
#' @import stats
#' @export
#' @examples
#' ### Number of groups of individuals to simulate ###
#' n_groups <- 2
#' ### Number of dates to simulate for each group ###
#' n_dates <- c(2, 3)
#' ### Setting up the parameters for the simulation ###
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of these delays
#' ### Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each group? ###
#' index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)))
#' ### Perform the simulation ###
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' ### Find the range ###
#' find_range(D$true_dat)
#' ### Compare with range specified in the first place ### 
#' range_dates
find_range <- function(obs_dat)
{
  min_date <- min(obs_dat[[1]][,1], na.rm=TRUE)
  max_date <- max(obs_dat[[1]][,1], na.rm=TRUE)
  for(g in 1:length(obs_dat))
  {
    for(j in 1:ncol(obs_dat[[g]]))
    {
      min_date_tmp <- min(obs_dat[[g]][,j], na.rm=TRUE)
      min_date <- min(c(min_date, min_date_tmp), na.rm=TRUE)
      
      max_date_tmp <- max(obs_dat[[g]][,j], na.rm=TRUE)
      max_date <- max(c(max_date, max_date_tmp), na.rm=TRUE)
    }
  }
  return(c(min_date, max_date))
}