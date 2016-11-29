###############################################
### move functions ###
###############################################

#' Performs one iteration of an MCMC move for the augmented data
#' 
#' @param i The index of the individual for whom augmented data should be moved
#' @param group_idx The index of the group for whom augmented data should be moved
#' @param date_idx The index of the date which should be moved
#' @param curr_aug_dat The current augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperpriors A list of hyperpriors: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}. If NULL, will be computed automatically.
#' @details \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperpriors} should be a list containing:
#' \itemize{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$CV}}
#' }
#' \code{index_dates} should be a list of length \code{n_groups=length(obs_dat)}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' The number of columns of index_dates[[k]] should match the length of theta$mu[[k]] and theta$CV[[k]] 
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are respectively the mean and coefficient of variation of two delays: the first delay being between date 1 and date 2, and the second being between date 1 and date 3. 
#' 
#' The function performs the move as follows, using a Metropolis algorithm. 
#' For the date to be moved, a new value is drawn from the marginal posterior of one of the delays this dates is involved in.
#' If the date is involved in several delays, one of the delays is randomly selected. 
#' The element E indicating whether the observed date is missing, recorded correctly or recorded with error, is adjusted accordingly given the proposed value of D. 
#' The new augmented data is then accepted with probability given by the ratio of the posterior values at the new augmented data and the old augmented data. 
#' @return A list of two elements:
#'  \itemize{
#'  \item{\code{new_aug_dat}}{: Same as \code{curr_aug_dat} but where the relevant dates have been updated}
#'  \item{\code{accept}}{: A scalar with value 1 if the move was accepted and 0 otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
move_Di <- function(i, group_idx, date_idx, 
                    curr_aug_dat,
                    theta, 
                    obs_dat, 
                    hyperpriors, 
                    index_dates,
                    range_dates=NULL) 
{
  if(is.null(range_dates)) range_dates <- find_range(obs_dat)
  
  # draw proposed value for D using one of the delays
  
  # which delays is this particular date involved in?
  
  x <- which(index_dates[[group_idx]]==date_idx, arr.ind = TRUE)
  which_delay <- x[2]
  from_idx <- sapply(1:nrow(x), function(k) index_dates[[group_idx]][-x[k,1],x[k,2]] )
  from_value <- sapply(1:nrow(x), function(k) curr_aug_dat$D[[group_idx]][i,index_dates[[group_idx]][-x[k,1],x[k,2]]])
  
  # if several delays involved, choose one at random
  tmp <- sample(1:length(from_idx), 1)
  from_idx <- from_idx[tmp]
  from_value <- from_value[tmp]
  param_delay <- find_params_gamma(theta$mu[[group_idx]][which_delay], CV=theta$CV[[group_idx]][which_delay])
  
  curr_aug_dat_value <- curr_aug_dat$D[[group_idx]][i,date_idx]
  sample_delay <- rgamma(1, shape=param_delay[1], scale=param_delay[2])
  if(date_idx<from_idx)
  {
    proposed_aug_dat_value <- from_value - round(sample_delay)
  }else
  {
    proposed_aug_dat_value <- from_value + round(sample_delay)
  }
  
  proposed_aug_dat <- curr_aug_dat
  proposed_aug_dat$D[[group_idx]][i,date_idx] <- proposed_aug_dat_value
  
  # adjust E_i accordingly
  # i.e. if D_i moves to y_i then E_i moves to 0, else E_i moves to 1. 
  missing <- which(is.na(obs_dat[[group_idx]][i,date_idx]))
  proposed_aug_dat$E[[group_idx]][i,date_idx][missing] <- -1 # y_i missing
  erroneous <- which(proposed_aug_dat$D[[group_idx]][i,date_idx]==obs_dat[[group_idx]][i,date_idx])
  proposed_aug_dat$E[[group_idx]][i,date_idx][erroneous] <- 0 # y_i observed without error
  non_erroneous <- which(!is.na(obs_dat[[group_idx]][i,date_idx]) & proposed_aug_dat$D[[group_idx]][i,date_idx]!=obs_dat[[group_idx]][i,date_idx] )
  proposed_aug_dat$E[[group_idx]][i,date_idx][non_erroneous] <- 1 # y_i observed with error
  
  # calculates probability of acceptance
  delay_idx <- which(index_dates[[group_idx]]==date_idx, arr.ind=TRUE)[,2] # these are the delays that are affected by the change in date date_idx
  ratio_post <- LL_observation_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates) - 
    LL_observation_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates) 
  ratio_post <- ratio_post + LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i) - 
    LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i)
  for(d in delay_idx)
    ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates) - 
    LL_delays_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates)
  
  ratio_post <- sum(ratio_post)
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, hyperpriors, index_dates) - 
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperpriors, index_dates)
  
  # no correction needed as this move is symetrical
  p_accept <- ratio_post 
  if(p_accept>0) {p_accept <- 0}
  
  # accept/reject step
  tmp <- log(runif(1))
  if(tmp<p_accept) # accepting with a certain probability
  {
    new_aug_dat <- proposed_aug_dat
    accept <- 1
  }else # reject
  {
    new_aug_dat <- curr_aug_dat
    accept <- 0
  }	
  
  # return a list of size 2 where 
  #		the first value is the new augmented data set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  return(list(new_aug_dat=new_aug_dat,accept=accept))
  
}
# test_move_Di <- move_Di(i=1, group_idx=1, date_idx=1, curr_aug_dat = aug_dat, theta, obs_dat, hyperpriors) 
# test_move_Di$new_aug_dat$D[[1]][1,1] # new value
# aug_dat$D[[1]][1,1] # old value

#' Performs one iteration of an MCMC move for either the parameter mu or the parameter CV (mean or CV of the various delays to be estimated)
#' 
#' @param what A string ("mu" or "CV") indicating which of the parameters to move
#' @param group_idx The index of the group for which mu or CV should be moved
#' @param delay_idx The index of the delay for which mu or CV should be moved
#' @param sdlog The standard deviation to be used for the proposal distribution, see details. 
#' @param aug_dat The augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param curr_theta The current list of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperpriors A list of hyperpriors: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @details \code{curr_theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperpriors} should be a list containing:
#' \itemize{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$CV}}
#' }
#' 
#' \code{index_dates} should be a list of length \code{n_groups=length(obs_dat)}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' The number of columns of index_dates[[k]] should match the length of theta$mu[[k]] and theta$CV[[k]] 
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are respectively the mean and coefficient of variation of two delays: the first delay being between date 1 and date 2, and the second being between date 1 and date 3. 
#' 
#' The function performs the move as follows, using a Metropolis-Hastings algorithm. 
#' For the parameter to be moved, a new value is drawn from a lognormal distribution with parameters \code{meanlog} equals the log of the current parameter value, and \code{sdlog=sdlog}.
#' The new parameter set is then accepted with probability given by the ratio of the posterior values at the new parameter set and the old parameter set, multiplie by a correction factor to reflect the non-symetrical nature of the move. 
#' @return A list of two elements:
#'  \itemize{
#'  \item{\code{new_theta}}{: Same as \code{curr_theta} but where \code{curr_theta$zeta} has been updated}
#'  \item{\code{accept}}{: A scalar with value 1 if the move was accepted and 0 otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
move_lognormal <- function(what=c("mu","CV"), group_idx, delay_idx, sdlog, 
                           aug_dat,
                           curr_theta, 
                           obs_dat, 
                           hyperpriors,
                           index_dates) 
{
  what <- match.arg(what)
  
  # draw proposed value
  curr_param_value <- curr_theta[[what]][[group_idx]][delay_idx]
  proposed_param_value <- rlnorm(1,meanlog=log(curr_param_value), sdlog=sdlog)
  
  proposed_theta <- curr_theta
  proposed_theta[[what]][[group_idx]][delay_idx] <- proposed_param_value
  
  # calculates probability of acceptance
  if(what=="mu")
  {
    ratio_post <- lprior_params_delay(what, proposed_theta, hyperpriors) - lprior_params_delay(what, curr_theta, hyperpriors) 
  }else if(what=="CV")
  {
    ratio_post <- lprior_params_delay(what, proposed_theta, hyperpriors) - lprior_params_delay(what, curr_theta, hyperpriors) 
  }
  Delta <- compute_delta_group_delay_and_indiv(aug_dat$D, group_idx, 1:nrow(obs_dat[[group_idx]]), delay_idx,  index_dates) # same for proposed and curent par values so no need to recompute twice
  ratio_post <- ratio_post + sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, proposed_theta, obs_dat, group_idx, delay_idx, 1:nrow(obs_dat[[group_idx]]), index_dates, Delta)) - 
    sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, curr_theta, obs_dat, group_idx, delay_idx, 1:nrow(obs_dat[[group_idx]]), index_dates, Delta)) 
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(aug_dat, proposed_theta, obs_dat, hyperpriors, index_dates) - 
  # lposterior_total(aug_dat, curr_theta, obs_dat, hyperpriors, index_dates)
  correction <- log(proposed_param_value) - log(curr_param_value) # correction for lognormal distribution
  p_accept <- ratio_post + correction # things are additive here as on log scale
  if(p_accept>0) {p_accept <- 0}
  
  # accept/reject step
  tmp <- log(runif(1))
  if(tmp<p_accept) # accepting with a certain probability
  {
    new_theta <- proposed_theta
    accept <- 1
  }else # reject
  {
    new_theta <- curr_theta 
    accept <- 0
  }	
  
  # return a list of size 2 where 
  #		the first value is the new parameter set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  return(list(new_theta=new_theta,accept=accept))
  
}
# test_move_mu <- move_lognormal(what="mu", group_idx=1, delay_idx=1, sdlog=0.1, aug_dat, curr_theta = theta, obs_dat, hyperpriors)
# test_move_mu$new_theta$mu[[1]][1] # new value
# theta$mu[[1]][1] # old value

#' Performs one iteration of an MCMC move for the parameter zeta (probability of a data being recorded erroneously, given it is recorded)
#' 
#' @param aug_dat The augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param curr_theta The current list of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperpriors A list of hyperpriors: see details.
#' @details \code{curr_theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperpriors} should be a list containing:
#' \itemize{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$CV}}
#' }
#' 
#' The function performs the move, using a Gibbs sampler. 
#' A new value of parameter zeta is drawn from its marginal posterior distribution, that is a beta distribution with parameters: 
#' \itemize{
#'  \item{\code{first shape parameter}}{: Equal to \code{hyperpriors$shape1_prob_error} + number_of_errors, where number_of_errors is the number of data points recorded with errors}
#'  \item{\code{second shape parameter}}{: Equal to \code{hyperpriors$shape2_prob_error} + number_of_recorded_dates-number_of_errors, where number_of_recorded_dates is the number of data points which are not missing, and number_of_errors is the number of data points recorded with errors}
#' }
#' @return A list of two elements:
#'  \itemize{
#'  \item{\code{new_theta}}{: Same as \code{curr_theta} but where \code{curr_theta$zeta} has been updated}
#'  \item{\code{accept}}{: A scalar with value 1 (as we are using a Gibbs sampler the move is always accepted)}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
move_zeta_gibbs <- function(aug_dat,
                            curr_theta, 
                            obs_dat, 
                            hyperpriors) 
{
  tmp <- compute_n_errors(aug_dat, obs_dat)
  number_of_errors <- tmp[1]
  number_of_recorded_dates <- tmp[2]
  
  # drawing from the marginal posterior distribution directly
  new_zeta <- rbeta(1, shape1=hyperpriors$shape1_prob_error + number_of_errors, shape2 = hyperpriors$shape2_prob_error + number_of_recorded_dates-number_of_errors)
  
  # therefore accept automatically 
  new_theta <- curr_theta
  new_theta$zeta <- new_zeta
  accept <- 1
  
  # return a list of size 2 where 
  #		the first value is the new parameter set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  return(list(new_theta=new_theta,accept=accept))
  
}
# test_move_zeta_gibbs <- move_zeta_gibbs(aug_dat, curr_theta = theta, obs_dat, hyperpriors)
# test_move_zeta_gibbs$new_theta$zeta # new value
# theta$zeta # old value
