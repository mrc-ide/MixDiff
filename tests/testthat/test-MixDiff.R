context("MixDiff") # one context by file, could be used to test different parts of the code, or either small unit test vs big behaviour things

test_that("Simplified likelihood ratio used in move_Di is correct", {
  
  set.seed(1)
  
  n_groups <- 4
  n_dates <- c(2, 3, 4, 4)
  
  theta <- list()
  theta$prop_missing_data <- 0.05 ### this is missing from the estimation model (directly available from the data so not explicitely modelled)
  theta$zeta <- 0.05 ### probability that, when not missing, the date is recorded with error
  theta$mu <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
  theta$CV <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))
  n_per_group <- rep(100, n_groups)
  range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
  index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
  index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )
  
  p_error <- list(external_swap=.04,internal_swap=.005,neighbour_substitution=0.05,distant_substitution=0.02,random=0.01)
  
  D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
  D_with_error <- simul_true_data(theta, n_per_group, range_dates, index_dates, simul_error = TRUE, p_error = p_error, remove_allNA_indiv=TRUE)
  tmp <- simul_obs_dat(D$true_dat, theta, range_dates)
  E <- tmp$E
  obs_dat <- tmp$obs_dat
  
  aug_dat <- list(D=D$true_dat, E=E)
  range_dates <- find_range(obs_dat) # update range based on observations
  
  # calculate the error matrix - for this we need actual dates not numbers
  rd <- int_to_date(range_dates)
  date_transition_mat_obs_true_log <- calculate_date_matrix(rd[1], rd[2], p_error, log = TRUE)
  
  hyperpriors <- list(
    shape1_prob_error=3, 
    shape2_prob_error=12, 
    mean_mean_delay=100, 
    mean_CV_delay=100)
  
  curr_aug_dat <- aug_dat
  
  i <- 1
  group_idx <- 1
  date_idx <- 1
  
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
  non_erroneous <- which(proposed_aug_dat$D[[group_idx]][i,date_idx]==obs_dat[[group_idx]][i,date_idx])
  proposed_aug_dat$E[[group_idx]][i,date_idx][non_erroneous] <- 0 # y_i observed without error
  erroneous <- which(!is.na(obs_dat[[group_idx]][i,date_idx]) & proposed_aug_dat$D[[group_idx]][i,date_idx]!=obs_dat[[group_idx]][i,date_idx] )
  proposed_aug_dat$E[[group_idx]][i,date_idx][erroneous] <- 1 # y_i observed with error
  
  # calculates probability of acceptance
  delay_idx <- which(index_dates[[group_idx]]==date_idx, arr.ind=TRUE)[,2] # these are the delays that are affected by the change in date date_idx
  ratio_post <- LL_observation_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates, date_transition_mat_obs_true_log = date_transition_mat_obs_true_log) - 
    LL_observation_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i, range_dates=range_dates, date_transition_mat_obs_true_log = date_transition_mat_obs_true_log) 
  ratio_post <- ratio_post + LL_error_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, date_idx, i) - 
    LL_error_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, date_idx, i)
  for(d in delay_idx)
    ratio_post <- ratio_post + LL_delays_term_by_group_delay_and_indiv(proposed_aug_dat, theta, obs_dat, group_idx, d, i, index_dates) - 
    LL_delays_term_by_group_delay_and_indiv(curr_aug_dat, theta, obs_dat, group_idx, d, i, index_dates)
  
  ratio_post <- sum(ratio_post)
  
  ### note that ratio_post should be the same as: 
  ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, hyperpriors, index_dates, range_dates=range_dates, date_transition_mat_obs_true_log = date_transition_mat_obs_true_log) - 
    lposterior_total(curr_aug_dat, theta, obs_dat, hyperpriors, index_dates, range_dates=range_dates, date_transition_mat_obs_true_log = date_transition_mat_obs_true_log)

})

test_that("Simplified likelihood ratio used in move_lognormal is correct", {
  
  set.seed(1)
  
  n_groups <- 4
  n_dates <- c(2, 3, 4, 4)
  
  theta <- list()
  theta$prop_missing_data <- 0.05 ### this is missing from the estimation model (directly available from the data so not explicitely modelled)
  theta$zeta <- 0.05 ### probability that, when not missing, the date is recorded with error
  theta$mu <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
  theta$CV <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))
  n_per_group <- rep(100, n_groups)
  range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
  index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
  index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )
  
  p_error <- list(external_swap=.04,internal_swap=.005,neighbour_substitution=0.05,distant_substitution=0.02,random=0.01)
  
  D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
  D_with_error <- simul_true_data(theta, n_per_group, range_dates, index_dates, simul_error = TRUE, p_error = p_error, remove_allNA_indiv=TRUE)
  tmp <- simul_obs_dat(D$true_dat, theta, range_dates)
  E <- tmp$E
  obs_dat <- tmp$obs_dat
  aug_dat <- list(D=D$true_dat, E=E)
  
  range_dates <- find_range(obs_dat) # update range based on observations
  
  # calculate the error matrix - for this we need actual dates not numbers
  rd <- int_to_date(range_dates)
  date_transition_mat_obs_true_log <- calculate_date_matrix(rd[1], rd[2], p_error, log = TRUE)
  
  hyperpriors <- list(
    shape1_prob_error=3, 
    shape2_prob_error=12, 
    mean_mean_delay=100, 
    mean_CV_delay=100)
  
  curr_theta <- theta
  
  for(what in c("mu", "CV"))
  {
    group_idx <- 1
    delay_idx <- 1
    sdlog <- 0.5
    
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
    ratio_post_long <- lposterior_total(aug_dat, proposed_theta, obs_dat, hyperpriors, index_dates, date_transition_mat_obs_true_log = date_transition_mat_obs_true_log) - 
      lposterior_total(aug_dat, curr_theta, obs_dat, hyperpriors, index_dates, date_transition_mat_obs_true_log = date_transition_mat_obs_true_log)
    
    expect_equal(ratio_post, ratio_post_long) # testing that the short version is equal to the full calculation
  }
})
