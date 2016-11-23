###############################################
### move functions ###
###############################################

## D_i ##
# probability of moving or not
# moves by drawing in the marginal posterior of one of the delays, and E_i is adjusted accordingly , 
# i.e. if D_i moves to y_i then E_i moves to 0, else E_i moves to 1. 
# then accept/reject based on posterior values
# this is symmetrical as independant so no correction needed
move_Di <- function(i, group_idx, date_idx, 
                    curr_aug_dat,
                    theta, 
                    obs_dat, 
                    shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_CV_delay=100, range_dates=NULL) 
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
  param_delay <- find_params_gamma_from_mean_CV(theta$mu[[group_idx]][which_delay], theta$CV[[group_idx]][which_delay])
  
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
    ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, d, i) - 
    LL_delays_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, d, i)
  
  ratio_post <- sum(ratio_post)
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay) - 
  # lposterior_total(curr_aug_dat, theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay)
  
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
# test_move_Di <- move_Di(i=1, group_idx=1, date_idx=1, curr_aug_dat = aug_dat, theta, obs_dat, shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_CV_delay=100) 
# test_move_Di$new_aug_dat$D[[1]][1,1] # new value
# aug_dat$D[[1]][1,1] # old value

### move mu with a lognormal proposal ### 
move_lognormal <- function(what=c("mu","CV"), group_idx, delay_idx, sdlog, 
                           aug_dat,
                           curr_theta, 
                           obs_dat, 
                           shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_CV_delay=100) 
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
    ratio_post <- lprior_params_delay(what, proposed_theta, prior_mean_mean_delay) - lprior_params_delay(what, curr_theta, prior_mean_mean_delay) 
  }else if(what=="CV")
  {
    ratio_post <- lprior_params_delay(what, proposed_theta, prior_mean_CV_delay) - lprior_params_delay(what, curr_theta, prior_mean_CV_delay) 
  }
  Delta <- compute_delta_group_delay_and_indiv(aug_dat$D, group_idx, delay_idx,  1:nrow(obs_dat[[group_idx]]), index = index_dates) # same for proposed and curent par values so no need to recompute twice
  ratio_post <- ratio_post + sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, proposed_theta, obs_dat, group_idx, delay_idx, 1:nrow(obs_dat[[group_idx]]), Delta)) - 
    sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, curr_theta, obs_dat, group_idx, delay_idx, 1:nrow(obs_dat[[group_idx]]), Delta)) 
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(aug_dat, proposed_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay) - 
  # lposterior_total(aug_dat, curr_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay)
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
# test_move_mu <- move_lognormal(what="mu", group_idx=1, delay_idx=1, sdlog=0.1, aug_dat, curr_theta = theta, obs_dat, shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_CV_delay=100)
# test_move_mu$new_theta$mu[[1]][1] # new value
# theta$mu[[1]][1] # old value

### move zeta with a truncated (<1) lognormal proposal ###
move_zeta_gibbs <- function(aug_dat,
                            curr_theta, 
                            obs_dat, 
                            shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_CV_delay=100) 
{
  tmp <- compute_n_errors(aug_dat, obs_dat)
  number_of_errors <- tmp[1]
  number_of_recorded_dates <- tmp[2]
  
  # drawing from the marginal posterior distribution directly
  new_zeta <- rbeta(1, shape1=shape1_prob_error + number_of_errors, shape2 = shape2_prob_error + number_of_recorded_dates-number_of_errors)
  
  # therefore accept automatically 
  new_theta <- curr_theta
  new_theta$zeta <- new_zeta
  accept <- 1
  
  # return a list of size 2 where 
  #		the first value is the new parameter set in the chain
  #		the second value is 1 if the proposed value was accepted, 0 otherwise
  return(list(new_theta=new_theta,accept=accept))
  
}
# test_move_zeta_gibbs <- move_zeta_gibbs(aug_dat, curr_theta = theta, obs_dat, shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_CV_delay=100)
# test_move_zeta_gibbs$new_theta$zeta # new value
# theta$zeta # old value
