###############################################
###############################################
### move functions ###
###############################################
###############################################

###############################################
### Move augmented dates D ###
###############################################

#' Performs one iteration of an MCMC move for the augmented data
#' 
#' @param i The index of the individual(s) for whom augmented data should be moved
#' @param group_idx The index of the group for whom augmented data should be moved
#' @param date_idx The index of the date which should be moved
#' @param curr_aug_dat The current augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}. If NULL, will be computed automatically.
#' @param tol A positive numerical value indicating the size of the tail of the distribution of different delays which can be ignored when drawing from these delays in the moves
#' @details \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
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
#' For the date to be moved, a new value is drawn from the marginal posterior of all the delays this date is involved in.
#' The element E indicating whether the observed date is missing, recorded correctly or recorded with error, is adjusted accordingly given the proposed value of D. 
#' The new augmented data is then accepted with probability given by the ratio of the posterior values at the new augmented data and the old augmented data, 
#' with a correction factor since the move is non symmetrical. 
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
                    hyperparameters, 
                    index_dates,
                    range_dates = NULL,
                    tol = 1e-6) 
{
  if(is.null(range_dates)) range_dates <- find_range(obs_dat)
  
  # draw proposed value for D using one of the delays
  # which delays is this particular date involved in?
  can_be_inferred_directly_from <- infer_directly_from(group_idx, date_idx, 
                                                       curr_aug_dat$D, i, theta) 
  
  # use this information to propose a clever D conditional on the current relevant delays parameters
  proposed <- sample_new_date_value(group_idx, 
                                    can_be_inferred_directly_from, 
                                    index_dates, tol = tol)
  tmp <- match(curr_aug_dat$D[[group_idx]][i, date_idx], 
               proposed$all_possible_values)
  if(!is.na(tmp))
  {
    proposed$probability_cur_value <- proposed$probabilities[tmp]
    
    proposed_aug_dat <- curr_aug_dat
    proposed_aug_dat$D[[group_idx]][i,date_idx] <- proposed$inferred
    
    # Now adjust E_i accordingly
    # i.e. if D_i moves to y_i then E_i moves to 0, else E_i moves to 1. 
    missing <- which(is.na(obs_dat[[group_idx]][i,date_idx]))
    proposed_aug_dat$E[[group_idx]][i,date_idx][missing] <- -1 # y_i missing
    non_erroneous <- which(proposed_aug_dat$D[[group_idx]][i,date_idx]==obs_dat[[group_idx]][i,date_idx])
    proposed_aug_dat$E[[group_idx]][i,date_idx][non_erroneous] <- 0 # y_i observed without error
    erroneous <- which(!is.na(obs_dat[[group_idx]][i,date_idx]) & proposed_aug_dat$D[[group_idx]][i,date_idx]!=obs_dat[[group_idx]][i,date_idx] )
    proposed_aug_dat$E[[group_idx]][i,date_idx][erroneous] <- 1 # y_i observed with error
    
    # calculates probability of acceptance
    delay_idx <- which(index_dates[[group_idx]]==date_idx, arr.ind=TRUE)[,2] # these are the delays that are affected by the change in date date_idx
    ratio_post <- LL_observation_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates) - 
      LL_observation_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates) 
    different_E <- proposed_aug_dat$E[[group_idx]][i,date_idx]!=curr_aug_dat$E[[group_idx]][i,date_idx]
    if(any(different_E)) # only need to look at the error term if some of the E have changed
    {
      ratio_post <- ratio_post + LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i) - 
        LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i)
    }
    for(d in delay_idx)
      ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates) - 
      LL_delays_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates)
    
    ratio_post <- sum(ratio_post)
    
    ### note that ratio_post should be the same as: 
    # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, hyperparameters, index_dates) - 
    #   lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates)
    # ratio_post_long <- sum(ratio_post_long)
    # if(abs(ratio_post - ratio_post_long)>1e-3)
    # {
    #   print("PROBLEM in move_Di ratio_post calculation\n")
    # }  
    
    # correction needed as this move is NOT symetrical
    prob_cur_to_new <- proposed$probability_inferred_value
    prob_new_to_cur <- proposed$probability_cur_value
    correction <- log(prob_new_to_cur) - log(prob_cur_to_new)
    #print(correction)
    
    p_accept <- ratio_post + correction
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
    
  }else # reject because the algorithm would never have proposed current value...
  {
    #print(paste("group", group_idx, "individual", i))
    new_aug_dat <- curr_aug_dat
    accept <- 0
  }
  
  # return a list of size 2 where 
  #		the first value is the new augmented data set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  return(list(new_aug_dat=new_aug_dat,accept=accept))
  
}
# test_move_Di <- move_Di(i=1, group_idx=1, date_idx=1, curr_aug_dat = aug_dat, theta, obs_dat, hyperparameters) 
# test_move_Di$new_aug_dat$D[[1]][1,1] # new value
# aug_dat$D[[1]][1,1] # old value

###############################################
### Move augmented indicator for whether date is correctly recorded E ###
###############################################

propose_move_from_E0_to_E1 <- function(i, group_idx, date_idx, 
                                       curr_aug_dat,
                                       theta, 
                                       obs_dat, 
                                       hyperparameters, 
                                       index_dates,
                                       range_dates=NULL, 
                                       tol = 1e-6)
{
  # which delays can this date be derived from
  can_be_inferred_directly_from <- infer_directly_from(group_idx, date_idx, 
                                                       curr_aug_dat$D, i, theta) 
  
  # use this information to propose a clever D conditional on the current relevant delays parameters
  proposed <- sample_new_date_value(group_idx, 
                                    can_be_inferred_directly_from, 
                                    index_dates, tol = tol)
  
  if(length(proposed$all_possible_values)<2 & proposed$inferred == obs_dat[[group_idx]][i,date_idx]) # then the only available date to draw from is the observed one and we can't propose this move from E0 to E1
  {
    # reject 
    return(NULL)
  }else
  {
    while(proposed$inferred == obs_dat[[group_idx]][i,date_idx]) ### whilst we haven't moved to a place where E=1, try again
    {
      proposed <- sample_new_date_value(group_idx, 
                                        can_be_inferred_directly_from, 
                                        index_dates, tol = tol)
    }
    proposed$probability_cur_value <- 1 # because going back to E0 there is only one choice
    
    return(proposed)
  }
}

compute_p_accept_move_from_E0_to_E1 <- function(i, group_idx, date_idx, 
                                                curr_aug_dat, proposed_aug_dat,
                                                prob_move = NULL,
                                                theta, 
                                                obs_dat, 
                                                hyperparameters, 
                                                index_dates,
                                                range_dates=NULL, 
                                                tol = 1e-6)
{
  proposed_aug_dat_value <- proposed_aug_dat$D[[group_idx]][i, date_idx]
  
  if(is.null(prob_move))
  {
    tmp <- propose_move_from_E0_to_E1(i, group_idx, date_idx, 
                                      curr_aug_dat,
                                      theta, 
                                      obs_dat, 
                                      hyperparameters, 
                                      index_dates,
                                      range_dates, 
                                      tol)
    if(any(tmp$all_possible_values == proposed_aug_dat_value))
    {
      prob_move <- tmp$probabilities[tmp$all_possible_values == proposed_aug_dat_value]
    } else
    {
      prob_move <- 0
    }
  }
  
  # calculates probability of acceptance
  delay_idx <- which(index_dates[[group_idx]] == date_idx, arr.ind=TRUE)[,2] # these are the delays that are affected by the change in date date_idx
  ratio_post <- LL_observation_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates) - 
    LL_observation_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates) 
  ratio_post <- ratio_post + LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i) - 
    LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i)
  for(d in delay_idx)
    ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates) - 
    LL_delays_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates)
  
  ratio_post <- sum(ratio_post)
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, hyperparameters, index_dates) - 
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates)
  # ratio_post_long <- sum(ratio_post_long)
  # if(ratio_post != ratio_post_long)
  # {
  #   print("PROBLEM")
  # }
  
  # this move is NOT symetrical --> correction needed
  prob_cur_to_new <- prob_move
  prob_new_to_cur <- 1 # going back to E0 there is only one choice
  #correction <- log(prob_new_to_cur) - log(prob_cur_to_new)
  correction <- - log(prob_cur_to_new)
  #print(correction)
  
  return(c(ratio_post, correction))
}

propose_move_from_E1_to_E0 <- function(i, group_idx, date_idx, 
                                       curr_aug_dat,
                                       theta, 
                                       obs_dat, 
                                       hyperparameters, 
                                       index_dates,
                                       range_dates)
{
  
  curr_aug_dat_value <- curr_aug_dat$D[[group_idx]][i,date_idx]
  
  proposed_aug_dat_value <- obs_dat[[group_idx]][i,date_idx]
  
  return(proposed_aug_dat_value)
}

compute_p_accept_move_from_E1_to_E0 <- function(i, group_idx, date_idx, 
                                                curr_aug_dat, proposed_aug_dat,
                                                theta, 
                                                obs_dat, 
                                                hyperparameters, 
                                                index_dates,
                                                range_dates, 
                                                tol = 1e-6)
{
  
  curr_aug_dat_value <- curr_aug_dat$D[[group_idx]][i, date_idx]
  
  # which delays can this date be derived from
  can_be_inferred_directly_from <- infer_directly_from(group_idx, date_idx, 
                                                       proposed_aug_dat$D, i, theta) 
  # use this information to propose a clever D conditional on the current relevant delays parameters
  # we only use this to check what would be the probability of the reverse move - to compute the correction factor correctly in MH
  proposed <- sample_new_date_value(group_idx, 
                                    can_be_inferred_directly_from, 
                                    index_dates, tol = tol)
  
  if(any(proposed$all_possible_values == curr_aug_dat_value))
  {
    prob_reverse_move <- proposed$probabilities[proposed$all_possible_values == curr_aug_dat_value]
  } else
  {
    prob_reverse_move <- 0
  }
  
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
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, hyperparameters, index_dates) - 
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates)
  
  prob_cur_to_new <- 1 # going to E0 there is only one choice
  prob_new_to_cur <- prob_reverse_move
  #correction <- log(prob_new_to_cur) - log(prob_cur_to_new)
  correction <- log(prob_new_to_cur)
  #print(correction)
  
  return(c( ratio_post,  correction))
}

#' Performs one iteration of an MCMC move for the augmented data representing the indicator of error in observations
#' 
#' @param i The index of the individual for whom augmented data should be moved
#' @param group_idx The index of the group for whom augmented data should be moved
#' @param date_idx The index of the date for which the error indicator which should be moved
#' @param curr_aug_dat The current augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}. If NULL, will be computed automatically.
#' @param tol A positive numerical value indicating the size of the tail of the distribution of different delays which can be ignored when drawing from these delays in the moves
#' @details \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
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
#' The function performs the move as follows, using a Metropolis-Hastings algorithm. 
#' If E=-1 nothing happens.
#' If E=1, we propose a move to E=0 and hence D=the observed data. 
#' If E=0, we propose a move to E=1. D is then moved as follows: a new value is drawn from the marginal posterior of all the delays this date is involved in, repeatdely until D falls on a different day than the observed date to be consistent with E=1. 
#' This move is not symmetrical so we use a correction factor in computing the probability of acceptance in the Metropolis Hastings which accounts for the asymetry. 
#' @return A list of two elements:
#'  \itemize{
#'  \item{\code{new_aug_dat}}{: Same as \code{curr_aug_dat} but where the relevant indicators of errors in dates have been updated}
#'  \item{\code{accept}}{: A scalar with value 1 if the move was accepted and 0 otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
move_Ei <- function(i, group_idx, date_idx, 
                    curr_aug_dat,
                    theta, 
                    obs_dat, 
                    hyperparameters, 
                    index_dates,
                    range_dates=NULL,
                    tol = 1e-6) 
{
  if(length(date_idx)>1)
  {
    print("potential issues move_Ei trying to move several dates at once")
  }
  if(length(i)>1) 
  {
    i <- i[1]
    warning("In move_Ei, i should be a numeric, not a vector. Using i[1] instead.")
  }
  
  if(is.null(range_dates)) range_dates <- find_range(obs_dat)
  
  curr_E_value <- curr_aug_dat$E[[group_idx]][i,date_idx]
  proposed_aug_dat <- curr_aug_dat
  
  if(curr_E_value!=-1) # if data not missing
  {
    new_E_value <- 1-curr_E_value
    proposed_aug_dat$E[[group_idx]][i,date_idx] <- new_E_value
    
    if(curr_E_value==0)# moving from E=0 to E=1
    {
      proposed <- propose_move_from_E0_to_E1(i, group_idx, date_idx, 
                                             curr_aug_dat,
                                             theta, 
                                             obs_dat, 
                                             hyperparameters, 
                                             index_dates,
                                             range_dates,
                                             tol = tol)
      
      if(!is.null(proposed))
      {
        proposed_aug_dat$D[[group_idx]][i,date_idx] <- proposed$inferred
        
        tmp <- compute_p_accept_move_from_E0_to_E1(i, group_idx, date_idx, 
                                                   curr_aug_dat, proposed_aug_dat,
                                                   prob_move = proposed$probability_inferred_value,
                                                   theta, 
                                                   obs_dat, 
                                                   hyperparameters, 
                                                   index_dates,
                                                   range_dates,
                                                   tol = tol)
        
        if(any(is.infinite(tmp))) p_accept <- -Inf else p_accept <- sum(tmp) 
        if(p_accept>0) p_accept <- 0
        
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
      } else # reject
      {
        new_aug_dat <- curr_aug_dat
        accept <- 0
      }
      
      # return a list of size 2 where 
      #		the first value is the new augmented data set in the chain
      #		the second value is 1 if the proposed value was accepted, 0 otherwise
      res <- list(new_aug_dat=new_aug_dat, accept=accept)
      
    }else if(curr_E_value==1)# moving from E=1 to E=0
    {
      proposed_aug_dat$D[[group_idx]][i,date_idx] <- propose_move_from_E1_to_E0(i, group_idx, date_idx, 
                                                                                curr_aug_dat,
                                                                                theta, 
                                                                                obs_dat, 
                                                                                hyperparameters, 
                                                                                index_dates,
                                                                                range_dates)
      
      tmp <- compute_p_accept_move_from_E1_to_E0(i, group_idx, date_idx, 
                                                 curr_aug_dat, proposed_aug_dat,
                                                 theta, 
                                                 obs_dat, 
                                                 hyperparameters, 
                                                 index_dates,
                                                 range_dates, tol = tol)
      if(any(is.infinite(tmp))) p_accept <- -Inf else p_accept <- sum(tmp)
      if(p_accept>0) p_accept <- 0
      
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
      res <- list(new_aug_dat=new_aug_dat, accept=accept)
    }
  }else # if E=-1, can't move. 
  {
    res <- list(new_aug_dat=curr_aug_dat, accept=0)
  }
  
  return(res)
}

###############################################
### Swap Es 
### (when only 2 Es, related by a delay, are recorded, 
### one with error and one without error, 
### propose to swap the two) ###
###############################################

find_2Eis_to_swap <- function(group_idx, curr_aug_dat)
{
  # find indexes of individuals in group_idx who, according to curr_aug_dat
  # have exactly one correctly recorded and one incorrectly recorded dates
  which(sapply(1:nrow(curr_aug_dat$E[[group_idx]]), function(i) 
    sum(curr_aug_dat$E[[group_idx]][i,]==1)==1 & 
      sum(curr_aug_dat$E[[group_idx]][i,]==0)==1))
}

find_Eis_to_swap <- function(group_idx, curr_aug_dat)
{
  # find indexes of individuals in group_idx who, according to curr_aug_dat
  # have at least one one correctly recorded and one incorrectly recorded dates
  which(sapply(1:nrow(curr_aug_dat$E[[group_idx]]), function(i) 
    length(unique(
      curr_aug_dat$E[[group_idx]][i,curr_aug_dat$E[[group_idx]][i,]!=-1]))>1))
}

#' Performs one iteration of an MCMC move for the augmented data where the indicators of error in observations for one individual are swapped, i.e. the errors become non errors and vice versa. 
#' 
#' @param i The index of the individual for whom augmented data should be moved
#' @param group_idx The index of the group for whom augmented data should be moved. This should be an individual and for whom in curr_aug_dat there is at least one correctly recorded date and one erroneously recorded dates.
#' @param curr_aug_dat The current augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param theta List of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @param range_dates A vector containing the range of dates in \code{obs_dat}. If NULL, will be computed automatically.
#' @param index_dates_order A list of same lenght as index_dates, containing indications on ordering of dates for each group - see \code{\link{compute_index_dates_order}}.
#' @param tol A positive numerical value indicating the size of the tail of the distribution of different delays which can be ignored when drawing from these delays in the moves
#' @details \code{theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
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
#' The function performs the move as follows, using a Metropolis-Hastings algorithm. 
#' It proposes to swap the indicators of errors in dates, e.g. 0 becomes 1 and 1 becomes 0.
#' First, for the dates where E = 1, we propose a move from E=1 to E=0, we automatically move D to the corresponding observed data. 
#' Then, for any dates not recorded (E = -1) we propose new values by drawing from the joint marginal posterior distribution of the relevant delays, using as reference dates only those that have just been moved to E=0 (not the ones which were initially at E=0 as these will move to erroneous values at the next step).
#' Finally, for the dates which were initially at E = 0, we propose a move from E=0 to E=1, D is then moved as follows: a new value is drawn from the joint marginal posterior of all the delays this date is involved in, repeatdely until D falls on a different day than the observed date to be consistent with E=1. 
#' If the date is involved in several delays, one of the delays is randomly selected. 
#' This move is not symmetrical so we use a correction factor in computing the probability of acceptance in the Metropolis Hastings which accounts for the asymetry. 
#' @return A list of two elements:
#'  \itemize{
#'  \item{\code{new_aug_dat}}{: Same as \code{curr_aug_dat} but where the relevant indicators of errors in dates have been updated}
#'  \item{\code{accept}}{: A scalar with value 1 if the move was accepted and 0 otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
swap_Ei <- function(i, group_idx,  
                    curr_aug_dat,
                    theta, 
                    obs_dat, 
                    hyperparameters, 
                    index_dates,
                    range_dates=NULL, 
                    index_dates_order = compute_index_dates_order(index_dates),
                    tol = 1e-6) 
{
  if(is.null(range_dates)) range_dates <- find_range(obs_dat)
  
  all_E_values <- curr_aug_dat$E[[group_idx]][i,]
  date_idx <- which(all_E_values %in% c(0,1))
  curr_E_values <- all_E_values[date_idx]
  
  date_idx_E0_to_E1 <- date_idx[curr_E_values==0] # currently correct date
  date_idx_E1_to_E0 <- date_idx[curr_E_values==1] # currently errouneous date
  
  # debug
  # if(group_idx==4 & i == 38)
  # {
  #  #browser()
  #  print("old state:")
  #  print(curr_aug_dat$E[[group_idx]][i,])
  #  print(curr_aug_dat$D[[group_idx]][i,])
  # }

  ###### First we move all E = 1 to E = 0 ######
  
  proposed_aug_dat_intermediate <- curr_aug_dat
  proposed_aug_dat_intermediate$E[[group_idx]][i,date_idx_E1_to_E0] <- 0
  proposed_aug_dat_intermediate$D[[group_idx]][i,date_idx_E1_to_E0] <- 
    propose_move_from_E1_to_E0(i, group_idx, date_idx_E1_to_E0, 
                               curr_aug_dat,
                               theta, 
                               obs_dat, 
                               hyperparameters, 
                               index_dates,
                               range_dates)
  
  # for the correction factor
  log_prob_move <- 0
  log_prob_reverse_move <- sum(sapply(1:length(date_idx_E1_to_E0), function(kk) 
  {
    tmp <- compute_p_accept_move_from_E0_to_E1(i, group_idx, date_idx_E1_to_E0[kk],
                                               proposed_aug_dat_intermediate, # the one we would be moving from 
                                               curr_aug_dat, # the one we would be moving to
                                               prob_move = NULL, # will be calculated inside function
                                               theta, obs_dat, hyperparameters, 
                                               index_dates,
                                               range_dates,
                                               tol = tol)
    return(-tmp[2])
  }))
  
   # if(group_idx==4 & i == 38)
   # {
   #   print("intermediate state 1:")
   #   print(proposed_aug_dat_intermediate$E[[group_idx]][i,])
   #   print(proposed_aug_dat_intermediate$D[[group_idx]][i,])
   #   #browser()
   # }
  
  ###### Then we update the missing values, if any, based on the data we just switched to 0 and ignoring the other ones (which will be switched to 1 later) ######
  
  ### The following code is to update D where E = -1 i.e. unrecorded dates
  ### linked to the date we've just changed
  missing_to_update <- which(
    proposed_aug_dat_intermediate$E[[group_idx]][i,] == -1)
  if(length(missing_to_update)>0)
  {
    #browser()
    tmp <- 
      infer_missing_dates(D = proposed_aug_dat_intermediate$D, 
                          E = proposed_aug_dat_intermediate$E, 
                          group_idx, i, index_dates_order, 
                          do_not_infer_from = date_idx_E0_to_E1, 
                          theta, tol = tol)
    proposed_aug_dat_intermediate$D[[group_idx]][i,missing_to_update] <- tmp$D[[group_idx]][i,missing_to_update] 
    prob <- tmp$prob
  }else
  {
    prob <- 1
  }
  
  # if(group_idx==4 & i == 38)
  # {
  #   print("intermediate state 2:")
  #   print(proposed_aug_dat_intermediate$E[[group_idx]][i,])
  #   print(proposed_aug_dat_intermediate$D[[group_idx]][i,])
  #   #browser()
  # }
  
  # for the correction factor ### TO DO: ADD THE CORRECT VALUES FOR CORRECTION FACTOR HERE
  log_prob_move <- c(log_prob_move, log(prob))
  # log_prob_reverse_move can not yet be calculated at this stage 
  # as it depends on where we would be moving from so need to know the final proposed_aug_dat
  
  ###### Finally we move the dates which were initially at E = 0 to E = 1 ######
  
  proposed_aug_dat <- proposed_aug_dat_intermediate
  proposed_aug_dat$E[[group_idx]][i,date_idx_E0_to_E1] <- 1
  
  for(kk in date_idx_E0_to_E1)
  {
    proposed <- propose_move_from_E0_to_E1(i, group_idx, kk, 
                                           proposed_aug_dat_intermediate,
                                           theta, 
                                           obs_dat, 
                                           hyperparameters, 
                                           index_dates,
                                           range_dates, 
                                           tol = tol)
    
    if(!is.null(proposed))
    {
      proposed_aug_dat$D[[group_idx]][i,kk] <- proposed$inferred
    } else
    {
      proposed_aug_dat$D[[group_idx]][i,kk] <- NULL
    }
  }
  
  # for the correction factor 
  if(any(is.null(proposed_aug_dat$D[[group_idx]][i,date_idx_E0_to_E1])))
  {
    log_prob_move <- log_prob_move(log_prob_move, -Inf)
  } else
  {
    log_prob_move <- c(log_prob_move, sum(sapply(1:length(date_idx_E0_to_E1), function(kk) 
    {
      tmp <- compute_p_accept_move_from_E1_to_E0(i, group_idx, date_idx_E0_to_E1[kk],
                                                 proposed_aug_dat, # the one we would be moving from 
                                                 proposed_aug_dat_intermediate, # the one we would be moving to
                                                 theta, obs_dat, hyperparameters, 
                                                 index_dates,
                                                 range_dates, 
                                                 tol = tol)
      return(tmp[2])
    })))
  }
  
  # if(group_idx==4 & i == 38)
  # {
  #   print("final proposed state:")
  #   print(proposed_aug_dat$E[[group_idx]][i,])
  #   print(proposed_aug_dat$D[[group_idx]][i,])
  #   browser()
  # }
  
  # probability of reverse move for missing data
  ### TO COMPUTE
  if(length(missing_to_update)>0)
  {
    #browser()
    tmp <- 
      infer_missing_dates(proposed_aug_dat_intermediate$D, 
                          E = proposed_aug_dat_intermediate$E, 
                          group_idx, i, index_dates_order, 
                          do_not_infer_from = date_idx_E1_to_E0, 
                          theta, 
                          tol = tol,
                          move_to_D = curr_aug_dat$D[[group_idx]][i,])
    prob_reverse <- tmp$prob_move_to_D
  }else
  {
    prob_reverse <- 1
  }
  log_prob_reverse_move <- c(log_prob_reverse_move, log(prob_reverse))
  # probability of reverse move for E0 to E1
  log_prob_reverse_move <- c(log_prob_reverse_move, 0) # as for the reverse move we go back to observed data E = 0 so only one choice
  
   # if(group_idx==4 & i == 38)
   # {
   #   print("final proposed state before calculation of p_accept:")
   #   print(proposed_aug_dat$E[[group_idx]][i,])
   #   print(proposed_aug_dat$D[[group_idx]][i,])
   #   browser()
   # }
  
  delay_idx <- which(colSums(matrix(index_dates[[group_idx]] %in% 
                                      date_idx_E1_to_E0, 
                                    nrow = nrow(index_dates[[group_idx]] )))>0)
  delay_idx <- c(delay_idx, 
                 which(colSums(matrix(index_dates[[group_idx]] %in% 
                                        date_idx_E0_to_E1, 
                                      nrow = 
                                        nrow(index_dates[[group_idx]] )))>0))
  delay_idx <- sort(unique(delay_idx))
  
  ratio_post_obs <- sum(LL_observation_term_by_group_delay_and_indiv(
    proposed_aug_dat, theta, obs_dat, group_idx, 
    date_idx_E1_to_E0, i, range_dates=range_dates) - 
      LL_observation_term_by_group_delay_and_indiv(
        curr_aug_dat, theta, obs_dat, group_idx, 
        date_idx_E1_to_E0, i, range_dates=range_dates) ) +  
    sum(LL_observation_term_by_group_delay_and_indiv(
      proposed_aug_dat, theta, obs_dat, group_idx, 
      date_idx_E0_to_E1, i, range_dates=range_dates) - 
        LL_observation_term_by_group_delay_and_indiv(
          curr_aug_dat, theta, obs_dat, group_idx, 
          date_idx_E0_to_E1, i, range_dates=range_dates) )
  
  #print(ratio_post_obs)
  ## should be the same as:
  # LL_observation_term(proposed_aug_dat, theta, obs_dat, range_dates) - LL_observation_term(curr_aug_dat, theta, obs_dat, range_dates)
  
  ratio_post_error <- sum(LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx_E1_to_E0, i) - 
                            LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx_E1_to_E0, i)) + 
    sum(LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx_E0_to_E1, i) - 
          LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx_E0_to_E1, i))
  #print(ratio_post_error)
  ## should be the same as:
  # LL_error_term(proposed_aug_dat, theta, obs_dat) - LL_error_term(curr_aug_dat, theta, obs_dat)
  ## LL_error_term_slow(proposed_aug_dat, theta, obs_dat) - LL_error_term_slow(curr_aug_dat, theta, obs_dat)
  
  ratio_post_delay <- 0
  for(d in delay_idx)
  {
    ratio_post_delay <- ratio_post_delay + sum(LL_delays_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates) - 
                                                 LL_delays_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates))
  }
  #print(ratio_post_delay)
  ## should be the same as:
  # LL_delays_term(proposed_aug_dat, theta, obs_dat, index_dates) - LL_delays_term(curr_aug_dat, theta, obs_dat, index_dates)
  
  ratio_post <- ratio_post_obs + ratio_post_error + ratio_post_delay
  
  ### should be the same as: 
  #ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, hyperparameters, index_dates) - 
  # lposterior_total(curr_aug_dat, theta, obs_dat, hyperparameters, index_dates)
  
  # This is not a symetric move so need a correction factor
  #print(paste("group", group_idx))
  #print(paste("individual", i))
  #if(group_idx == 4 & i ==21)
  #{
  #  browser()
  #}
  # if(length(missing_to_update)>0)
  # {
  #   browser()
  # }
  p_accept <- ratio_post + sum(log_prob_reverse_move) - sum(log_prob_move)
  #print(p_accept)
  if(p_accept>0) p_accept <- 0
  
  # if(group_idx==4 & i == 38)
  # {
  #   browser()
  # }
  
  # accept/reject step
  tmp <- log(runif(1))
  if(tmp<p_accept) # accepting with a certain probability
  {
    new_aug_dat <- proposed_aug_dat
    accept <- 1
    # if(group_idx == 4 && i == 5)
    # {
    #   print("accepted")
    # }
    # if(length(missing_to_update)>0)
    # {
    #   print("accepted swap with missing values")
    # }
  }else # reject
  {
    new_aug_dat <- curr_aug_dat
    accept <- 0
    # if(group_idx == 4 && i == 5)
    # {
    #   print("rejected")
    # }
  }	
  
  # return a list of size 2 where 
  #		the first value is the new augmented data set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  res <- list(new_aug_dat=new_aug_dat, accept=accept)
  
  return(res)
  
}

###############################################
### Move mean or CV of delay ###
###############################################

#' Performs one iteration of an MCMC move for either the parameter mu or the parameter CV (mean or CV of the various delays to be estimated)
#' 
#' @param what A string ("mu" or "CV") indicating which of the parameters to move
#' @param group_idx The index of the group for which mu or CV should be moved
#' @param delay_idx The index of the delay for which mu or CV should be moved
#' @param sdlog The standard deviation to be used for the proposal distribution, see details. 
#' @param aug_dat The augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param curr_theta The current list of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @details \code{curr_theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
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
                           hyperparameters,
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
    ratio_post <- lprior_params_delay(what, proposed_theta, hyperparameters) - lprior_params_delay(what, curr_theta, hyperparameters) 
  }else if(what=="CV")
  {
    ratio_post <- lprior_params_delay(what, proposed_theta, hyperparameters) - lprior_params_delay(what, curr_theta, hyperparameters) 
  }
  Delta <- compute_delta_group_delay_and_indiv(aug_dat$D, group_idx, seq_len(nrow(obs_dat[[group_idx]])), delay_idx,  index_dates) # same for proposed and curent par values so no need to recompute twice
  ratio_post <- ratio_post + sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, proposed_theta, obs_dat, group_idx, delay_idx, seq_len(nrow(obs_dat[[group_idx]])), index_dates, Delta)) - 
    sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, curr_theta, obs_dat, group_idx, delay_idx, seq_len(nrow(obs_dat[[group_idx]])), index_dates, Delta)) 
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(aug_dat, proposed_theta, obs_dat, hyperparameters, index_dates) - 
  # lposterior_total(aug_dat, curr_theta, obs_dat, hyperparameters, index_dates)
  
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
# test_move_mu <- move_lognormal(what="mu", group_idx=1, delay_idx=1, sdlog=0.1, aug_dat, curr_theta = theta, obs_dat, hyperparameters)
# test_move_mu$new_theta$mu[[1]][1] # new value
# theta$mu[[1]][1] # old value

###############################################
### Move zeta (probability of erroneous recording of dates) ###
###############################################

#' Performs one iteration of an MCMC move for the parameter zeta (probability of a data being recorded erroneously, given it is recorded)
#' 
#' @param aug_dat The augmented data; a list of observed data, in the format returned by \code{\link{simul_true_data}}. 
#' @param curr_theta The current list of parameters; see details.
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param hyperparameters A list of hyperparameters: see details.
#' @details \code{curr_theta} should be a list containing:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups} (the number of groups to be simulated data). Each element of \code{mu} should be a scalar of vector giving the mean delay(s) to use for simulation of dates in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar of vector giving the coefficient o variation of the delay(s) to use for simulation of dates in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' \code{hyperparameters} should be a list containing:
#' \itemize{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$CV}}
#' }
#' 
#' The function performs the move, using a Gibbs sampler, which we can use because 
#' we use a bernoulli likelihood for the presence/absence of error for each date, 
#' and we chose a conjugate beta prior for the probability of erroneous entry. 
#' A new value of parameter zeta is drawn from its marginal posterior distribution, 
#' which is therefore a beta distribution with parameters: 
#' \itemize{
#'  \item{\code{first shape parameter}}{: Equal to \code{hyperparameters$shape1_prob_error} + number_of_errors, where number_of_errors is the number of data points recorded with errors}
#'  \item{\code{second shape parameter}}{: Equal to \code{hyperparameters$shape2_prob_error} + number_of_recorded_dates-number_of_errors, where number_of_recorded_dates is the number of data points which are not missing, and number_of_errors is the number of data points recorded with errors}
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
                            hyperparameters) 
{
  tmp <- compute_n_errors(aug_dat, obs_dat)
  number_of_errors <- tmp[1]
  number_of_recorded_dates <- tmp[2]
  
  # drawing from the marginal posterior distribution directly
  new_zeta <- rbeta(1, 
                    shape1 = hyperparameters$shape1_prob_error + 
                      number_of_errors, 
                    shape2 = hyperparameters$shape2_prob_error + 
                      number_of_recorded_dates - number_of_errors)
  
  # therefore accept automatically 
  new_theta <- curr_theta
  new_theta$zeta <- new_zeta
  accept <- 1
  
  # return a list of size 2 where 
  #		the first value is the new parameter set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  return(list(new_theta=new_theta,accept=accept))
  
}
# test_move_zeta_gibbs <- move_zeta_gibbs(aug_dat, curr_theta = theta, obs_dat, hyperparameters)
# test_move_zeta_gibbs$new_theta$zeta # new value
# theta$zeta # old value
