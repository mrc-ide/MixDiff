#######################################
### functions to simulate a dataset ###
#######################################

#' Simulates data; see details
#' 
#' @param theta List of parameters; see details.
#' @param n_per_group Vector containing the number of individuals to simulate in each group
#' @param range_dates Range of integers in which to draw the first set of dates (these will ve drawn unifromly in that range)
#' @param index_dates A list containing indications on which delays to consider in the simulation, see details.
#' @param simul_error A boolean indicating whether to also simulate missingness and error in data or not (also see \code{\link[MixDiff]{simul_obs_dat}}).
#' @param remove_allNA_indiv A boolean stating whether individuals with only missing observations should be removed or not; only used if \code{simul_error} is TRUE (also see \code{\link[MixDiff]{simul_obs_dat}}).
#' @details \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{prop_missing_data} (only required if \code{simul_error} is TRUE)}{: A scalar in [0;1] giving the probability of each data point being missing.}
#'  \item{\code{zeta} (only required if \code{simul_error} is TRUE)}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{n_per_group} should be a vector of length \code{n_groups}.
#' 
#' \code{index_dates} should be a list of length \code{n_groups}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' The number of columns of index_dates[[k]] should match the length of theta$mu[[k]] and theta$CV[[k]] 
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are respectively the mean and coefficient of variation of two delays: the first delay being between date 1 and date 2, and the second being between date 1 and date 3. 
#' In the simulation, date 1 will be drawn uniformly within \code{range_dates}. 
#' Then date 2 will be drawn as date 1 + a discretised gamma distribution with mean theta$mu[[k]][1] and theta$CV[[k]][1]. 
#' Finally, date 3 will be drawn as date 1 + a discretised gamma distribution with mean theta$mu[[k]][1] and theta$CV[[k]][1]. 
#' @return A list of three items. 
#'  \itemize{
#'  \item{\code{true_dat}}{ A list of length \code{length(n_per_group)} matrices; each has \code{length(n_per_group)} rows corresponding to individuals and a certain number of columns derived from \code{index_dates}. 
#'  Elements of the matrices are integers corresponding to dates (see \code{\link[MixDiff]{int_to_date}} and \code{\link[MixDiff]{date_to_int}})}
#'  \item{\code{obs_dat}}{ A list structured as \code{true_dat} but where missing data and errors have been introduced (NULL if \code{simul_error} is FALSE)}
#'  \item{\code{E}} { A list structured similarly to \code{true_dat} and \code{obs_dat}, containing indicators of where \code{obs_dat} is missing (\code{E=-1}), 
#'  where \code{obs_dat} is recorded but with error (\code{E=1}), and where \code{obs_dat} is recorded with no error (\code{E=0}) (NULL if \code{simul_error} is FALSE)}
#'  }
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
simul_true_data <- function(theta, n_per_group, range_dates, index_dates, simul_error=FALSE, remove_allNA_indiv=FALSE)
{
  D <- list() 
  for(g in seq_len(length(theta$mu)) )
  {
    D[[g]] <- matrix(NA, n_per_group[g], length(theta$mu[[g]])+1)
    D[[g]][,1] <- sample(seq(range_dates[1],range_dates[2],1), n_per_group[g], replace = TRUE)
    for(j in seq_len(ncol(index_dates[[g]])) )
    {
      params <- find_params_gamma(theta$mu[[g]][j], CV=theta$CV[[g]][j])
      delay <- rgamma(n_per_group[g], shape=params[1], scale=params[2])
      D[[g]][,index_dates[[g]][2,j]]  <- D[[g]][,index_dates[[g]][1,j]] + round(delay)
    }
  }
  if(simul_error)
  {
    observed_D <- simul_obs_dat(D, theta, range_dates)
    return(list(true_dat=D, obs_dat=observed_D$obs_dat, E=observed_D$E))
  }else{
    return(list(true_dat=D, obs_dat=NULL, E=NULL))
  }
}

#' Introduces missingness and errors in data
#' 
#' @param D A list of data, in the format of the first element (called \code{true_dat}) in the list returned by \code{\link{simul_true_data}}. 
#' @param theta A list of parameters; see details.
#' @param range_dates Range of integers in which to draw the erroneous data (these will ve drawn unifromly in that range)
#' @param remove_allNA_indiv A boolean stating whether individuals with only NA dates should be removed or not. 
#' @details \code{theta} should be a list containing
#' \itemize{
#'  \item{\code{prop_missing_data}}{: A scalar in [0;1] giving the probability of each data point being missing.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' @return A list with two elements: 
#' \itemize{
#'  \item{\code{obs_dat}}{: A list similar to \code{D}, but where some data points are now missing, and some are erroneous}
#'  \item{\code{E}}{: A list structured similarly to \code{D} and \code{obs_dat}, containing indicators of where \code{obs_dat} is missing (\code{E=-1}), where \code{obs_dat} is recorded but with error (\code{E=1}), and where \code{obs_dat} is recorded with no error (\code{E=0})}
#' }
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
#' theta$prop_missing_data <- 0.25 # probability of data missing in observations
#' theta$zeta <- 0.05 # probability that, when not missing, the date is recorded with error
#' ### Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each group? ###
#' index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)))
#' ### Perform the simulation ###
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' observed_D <- simul_obs_dat(D$true_dat, theta, range_dates, remove_allNA_indiv=TRUE)
#' ### the observed dataset is smaller than the true one 
#' ### (some individuals with only missing data do not appear in the observed dataset)
#' nrow(D$true_dat[[1]])
#' nrow(observed_D$obs_dat[[1]])
simul_obs_dat <- function(D, theta, range_dates, remove_allNA_indiv=FALSE)
{
  E <- D
  obs_dat <- D
  for(g in seq_len(length(D)) )
  {
    for(j in seq_len(ncol(D[[g]])) )
    {
      E[[g]][,j] <- sample(c(-1, 1, 0), nrow(D[[g]]), replace=TRUE, prob=c(theta$prop_missing_data, (1-theta$prop_missing_data)*theta$zeta, (1-theta$prop_missing_data)*(1-theta$zeta)))
      obs_dat[[g]][E[[g]][,j]==-1,j]  <- NA
      obs_dat[[g]][E[[g]][,j]==0,j]  <- D[[g]][E[[g]][,j]==0,j]
      obs_dat[[g]][E[[g]][,j]==1,j]  <- sample(seq(range_dates[1], range_dates[2], 1), sum(E[[g]][,j]==1), replace = TRUE) # need to update if change error model
    }
    if(remove_allNA_indiv)
    {
      exclude <- which(rowSums(is.na(obs_dat[[g]]))==ncol(obs_dat[[g]]))
      if(length(exclude)>0)
      {
        obs_dat[[g]] <- obs_dat[[g]][-exclude,]
        E[[g]] <- E[[g]][-exclude,]
      }
    }
  }
  return(list(obs_dat=obs_dat, E=E))
}


