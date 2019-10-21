###############################################
### define augmented data to be used for initialisation of the chain ###
###############################################

are_dates_incompatible <- function(date1, date2, mindelay, maxdelay)
{
  return (date2-date1<mindelay | date2-date1>maxdelay)
}

### D contains the unobserved true dates ###

#' Initialises augmented data based on observed data (for the MCMC)
#' 
#' @param obs_dat A list of data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @param MCMC_settings A list of settings to be used for initialising the augmented data for the MCMC, see details.
#' @details \code{index_dates} should be a list; each elements corresponding to a group of individuals of interest. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and c(1, 3), this indicates that for group \code{k} we are interested in two delays: the first delay being between date 1 and date 2, and the second being between date 1 and date 3. 
#' \code{MCMC_settings} should be a list containing:
#' \itemize{
#'  \item{\code{init_options}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{mindelay}}{: The minimum delay, below which dates are considered incompatile with one another at the initialisation stage of the MCMC.}
#'  \item{\code{maxdelay}}{: The maximum delay, above which dates are considered incompatile with one another at the initialisation stage of the MCMC.  }
#'  \item{\code{record_every}}{: A number indicating, after the burnin, every how many iterations outputs should be recorded.}
#'  }
#'  }
#' }
#' @return A list with two elements: 
#' \itemize{
#'  \item{\code{D}}{: A list similar to \code{obs_dat}, but where no data points are missing, and some dates have been corrected to be consistent with the ordering rules inherent to \code{index_dates}}
#'  \item{\code{E}}{: A list structured similarly to \code{D} and \code{obs_dat}, containing indicators of where \code{obs_dat} is missing (\code{E=-1}), where \code{obs_dat} is recorded but with error (\code{E=1}), and where \code{obs_dat} is recorded with no error (\code{E=0})}
#' }
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
#' ### Initialise augmented data ###
#' MCMC_settings <- list(init_options=list(mindelay=0, maxdelay=100))
#' aug_dat <- initialise_aug_data(observed_D$obs_dat, index_dates, MCMC_settings)
initialise_aug_data <- function(obs_dat, index_dates, MCMC_settings)
{
  index_dates_order <- compute_index_dates_order(index_dates)
  n_groups <- length(obs_dat)
  D <- list()
  for(g in seq_len(n_groups) )
  {
    #print(paste("group", g))
    D[[g]] <- obs_dat[[g]]
    for(e in seq_len(nrow(D[[g]])))
    {
      #print(paste("individual", e))
      #print(D[[g]][e,])
      # first deal with incompatible dates
      for(j in seq_len(ncol(index_dates_order[[g]])))
      {
        #print(paste("date", j))
        if(!any(is.na(D[[g]][e,index_dates_order[[g]][,j]])))
        {
          # there is a problem if the dates have too short or too long delay
          #print(D[[g]][e,])
          if(are_dates_incompatible(D[[g]][e,index_dates_order[[g]][1,j]], D[[g]][e,index_dates_order[[g]][2,j]], MCMC_settings$init_options$mindelay, MCMC_settings$init_options$maxdelay) )
          {
            # check if there is one of the dates involved in more than one problematic delays, if so must be the problematic one:
            tmp <- table(as.vector(index_dates_order[[g]][,sapply(seq_len(ncol(index_dates_order[[g]])), function(j) are_dates_incompatible(D[[g]][e,index_dates_order[[g]][1,j]], D[[g]][e,index_dates_order[[g]][2,j]], MCMC_settings$init_options$mindelay, MCMC_settings$init_options$maxdelay) )]))
            if(any(tmp>1))
            {
              must_be_wrong <- which.max(tmp)[1]
            }else
            {
              # check which of all dates is most outlier compared to all other dates, and if several take the first one as the wrong one
              diff_from_median <- abs(D[[g]][e,index_dates_order[[g]][,j]] - median(D[[g]][e,], na.rm=TRUE))
              must_be_wrong <- which(diff_from_median %in% max(diff_from_median))[1]
              must_be_wrong <- index_dates_order[[g]][,j][must_be_wrong]
            }
            D[[g]][e,must_be_wrong] <- NA
            while(!(must_be_wrong %in% index_dates_order[[g]][,j]))
            {
              # check if there is one of the dates involved in more than one problematic delays, if so must be the problematic one:
              tmp <- table(as.vector(index_dates_order[[g]][,sapply(seq_len(ncol(index_dates_order[[g]])), function(j) are_dates_incompatible(D[[g]][e,index_dates_order[[g]][1,j]], D[[g]][e,index_dates_order[[g]][2,j]], MCMC_settings$init_options$mindelay, MCMC_settings$init_options$maxdelay) )]))
              if(any(tmp>1))
              {
                must_be_wrong <- which.max(tmp)[1]
              }else
              {
                # check which of all dates is most outlier compared to all other dates, and if several take the first one as the wrong one
                diff_from_median <- abs(D[[g]][e,index_dates_order[[g]][,j]] - median(D[[g]][e,], na.rm=TRUE))
                must_be_wrong <- which(diff_from_median %in% max(diff_from_median))[1]
                must_be_wrong <- index_dates_order[[g]][,j][must_be_wrong]
              }
              D[[g]][e,must_be_wrong] <- NA
            }
          }
        }
      }
      
      # now deal with missing dates
      D <- infer_missing_dates(D, E = NULL, g, e, index_dates_order)$D
      
    }
  }
  names(D) <- names(obs_dat)
  
  # compute E accordingly
  E <- list()
  for(g in seq_len(n_groups) )
  {
    E[[g]] <- matrix(NA,nrow(obs_dat[[g]]),ncol(obs_dat[[g]]))
    for(j in seq_len(ncol(obs_dat[[g]])) )
    {
      error <- D[[g]][,j] != obs_dat[[g]][,j]
      E[[g]][which(error),j] <- 1 # error
      E[[g]][which(!error),j] <- 0 # no error
      E[[g]][which(is.na(obs_dat[[g]][,j])), j] <- -1 # missing value
    }
    names(E[[g]]) <- names(obs_dat[[g]])
  }
  names(E) <- names(obs_dat)
  
  aug_dat <- list(D = D,
                  E = E)
  
  return(aug_dat)
}

###############################################
### define parameters to be used for initialisation of the chain ###
###############################################

#' Initialises parameters data on augmented data (for the MCMC)
#' 
#' @param aug_dat A list of data, in the format returned by \code{\link{simul_true_data}}. 
#' @param index_dates A list containing indications on which delays to consider in the simulation, same as in \code{\link{simul_true_data}}.
#' @param zeta_init A scalar giving the value zeta should be initialised to. 
#' @return A list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
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
#' ### Initialise augmented data first ###
#' MCMC_settings <- list(init_options=list(mindelay=0, maxdelay=100))
#' aug_dat <- initialise_aug_data(observed_D$obs_dat, index_dates, MCMC_settings)
#' ### Now initialise parameters based on the augmented data above ###
#' theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)
initialise_theta_from_aug_dat <- function(aug_dat, index_dates, zeta_init=0.1) # zeta_init doesn't really matter as we then use Gibbs sampler so will move fast to better values
{
  n_groups <- length(aug_dat$D)
  n_dates <- sapply(aug_dat$D, ncol)
    
  ### mean and std of distribution of various delays, by group
  ### we use a the starting point the observed mean and std of each delay in each group
  obs_delta <- compute_delta(aug_dat$D, index_dates)
  mu <- lapply(seq_len(n_groups), function(g) abs(apply(obs_delta[[g]], 2, mean, na.rm=TRUE) ))
  names(mu) <- names(n_dates)
  sigma <- lapply(seq_len(n_groups), function(g) abs(apply(obs_delta[[g]], 2, sd, na.rm=TRUE) ))
  names(sigma) <- names(n_dates)
  CV <- lapply(seq_len(n_groups), function(g) sigma[[g]]/mu[[g]])
  names(CV) <- names(n_dates)
  
  ### list of all parameters
  theta <- list(zeta = zeta_init, # zeta is the probability for a date to be misrecorded, conditional on being recorded (<-> Ei != - 1)
                # TODO:
                # could consider having zeta being type of date specific (e.g. more error on onset than death dates),
                # time specific and/or space specific
                mu = mu, # mean of gamma distributions used to characterise the various delays in different groups: mu[[g]][k] is the mean k^th delay in group g
                CV = CV) # CV of gamma distributions used to characterise the various delays in different groups: CV[[g]][k] is the CV k^th delay in group g
  
  return(theta)
  
}
