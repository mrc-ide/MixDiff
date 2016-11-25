#######################################
### functions to simulate a dataset ###
#######################################

#' Simulates data; see details
#' 
#' @param theta list of parameters; see details.
#' @param n_per_group Vector containing the number of individuals to simulate in each group
#' @param range_dates Range of integers in which to draw the first set of dates (these will ve drawn unifromly in that range)
#' @param index_dates A list containing indications on which delays to consider in the simulation, see details.
#' @details \code{theta} should be a list containing
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
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
#' @return A list of length \code{length(n_per_group)} matrices; each has \code{length(n_per_group)} rows corresponding to individuals and a certain number of columns derived from \code{index_dates}. Elements of the matrices are integers corresponding to dates (see \code{\link[MixDiff]{int_to_date}} and \code{\link[MixDiff]{date_to_int}})
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
simul_true_data <- function(theta, n_per_group, range_dates, index_dates)
{
  D <- list() 
  for(g in 1:length(theta$mu))
  {
    D[[g]] <- matrix(NA, n_per_group[g], length(theta$mu[[g]])+1)
    D[[g]][,1] <- sample(range_dates[1]:range_dates[2], n_per_group[g], replace = TRUE)
    for(j in 1:ncol(index_dates[[g]]))
    {
      params <- find_params_gamma(theta$mu[[g]][j], CV=theta$CV[[g]][j])
      delay <- rgamma(n_per_group[g], shape=params[1], scale=params[2])
      D[[g]][,index_dates[[g]][2,j]]  <- D[[g]][,index_dates[[g]][1,j]] + round(delay)
    }
  }
  return(D)
}

simul_obs_dat <- function(D, theta, range_dates)
{
  E <- D
  obs_dat <- D
  for(g in 1:length(D))
  {
    for(j in 1:ncol(D[[g]]))
    {
      E[[g]][,j] <- sample(c(-1, 1, 0), nrow(D[[g]]), replace=TRUE, prob=c(theta$prop_missing_data, (1-theta$prop_missing_data)*theta$zeta, (1-theta$prop_missing_data)*(1-theta$zeta)))
      obs_dat[[g]][E[[g]][,j]==-1,j]  <- NA
      obs_dat[[g]][E[[g]][,j]==0,j]  <- D[[g]][E[[g]][,j]==0,j]
      obs_dat[[g]][E[[g]][,j]==1,j]  <- sample(range_dates[1]:range_dates[2], sum(E[[g]][,j]==1), replace = TRUE) # need to update if change error model
    }
  }
  return(list(obs_dat=obs_dat, E=E))
}


