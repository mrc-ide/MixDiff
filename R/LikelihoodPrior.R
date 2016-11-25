###############################################
### likelihood function ###
###############################################

LL_observation_term_by_group_delay_and_indiv <- function(aug_dat, theta, obs_dat, group_idx, date_idx, indiv_idx, range_dates=NULL)
{
  if(is.null(range_dates)) range_dates <- find_range(obs_dat)
  LL <- vector()
  ### making sure D=y if E=0 ### note could remove this if by construction this is always true - could speed up code
  no_error <- which(aug_dat$E[[group_idx]][indiv_idx, date_idx] %in% 0)
  LL[no_error] <- log(aug_dat$D[[group_idx]][indiv_idx, date_idx][no_error] == obs_dat[[group_idx]][indiv_idx, date_idx][no_error]) 
  ### if E=1, what is the relationship between true date D and observed date y
  # for now, observation likelihood conditional on E=1 is uniform on the range of observed dates
  ### same for E=-1, D can take any value in range with same probability as they are all consistent with y=NA
  error_or_missing <- which((aug_dat$E[[group_idx]][indiv_idx, date_idx] %in% -1 | aug_dat$E[[group_idx]][indiv_idx, date_idx] %in% 1))
  ### K is the relative probability of observing a given error, conditional on presence of error
  # for now, K is given as 1/n, where n is the number of dates in the range_dates
  # could use something different if we define the space of possible errors differently. 
  ### think about impact of this choice
  K <- (1/as.numeric(diff(range_dates)))
  LL[error_or_missing] <- log( K* ((aug_dat$D[[group_idx]][indiv_idx, date_idx][error_or_missing] >= range_dates[1]) & (aug_dat$D[[group_idx]][indiv_idx, date_idx][error_or_missing] <= range_dates[2])) ) 
  
  LL[is.infinite(LL)] <- -100000 # arbitrarily small number to avoid -Inf
  
  return(LL)
}

LL_observation_term<-function(aug_dat, theta, obs_dat, range_dates=NULL)
{
  if(is.null(range_dates)) range_dates <- find_range(obs_dat)
  LL <- sum (sapply(1:length(obs_dat), function(g) sum (sapply(1:ncol(aug_dat$D[[g]]), function(j) sum(LL_observation_term_by_group_delay_and_indiv(aug_dat, theta, obs_dat, g, j, 1:nrow(obs_dat[[g]]),range_dates)) ) ) ) )
  return(LL)
}
# LL_observation_term(aug_dat, theta, obs_dat)

LL_error_term_by_group_delay_and_indiv <- function(aug_dat, theta, obs_dat, group_idx, date_idx, indiv_idx)
{
  res <- vector()
  missing <- which(aug_dat$E[[group_idx]][indiv_idx,date_idx] %in% -1)
  non_missing <- which(!(aug_dat$E[[group_idx]][indiv_idx,date_idx]  %in% -1))
  res[non_missing] <- log(theta$zeta)*as.numeric(aug_dat$E[[group_idx]][indiv_idx[non_missing],date_idx]) + log(1-theta$zeta)*(1-as.numeric(aug_dat$E[[group_idx]][indiv_idx[non_missing],date_idx]))
  res[missing] <- 0
  return(res)
}

compute_n_errors <- function(aug_dat, obs_dat)
{
  n_groups <- length(obs_dat)
  number_of_errors <- sum(sapply(1:n_groups, function(g) sum(aug_dat$E[[g]] %in% 1)))
  number_of_recorded_dates <- sum(sapply(1:n_groups, function(g) sum(!(aug_dat$E[[g]] %in% -1))))
  return(c(number_of_errors, number_of_recorded_dates))
}
# system.time(compute_n_errors(aug_dat, obs_dat))

LL_error_term<-function(aug_dat, theta, obs_dat)
{
  tmp <- compute_n_errors(aug_dat, obs_dat)
  number_of_errors <- tmp[1]
  number_of_recorded_dates <- tmp[2]
  #result<-dbinom(number_of_errors,number_of_recorded_dates,theta$zeta,log=TRUE)
  result<- log(theta$zeta)*number_of_errors + log(1-theta$zeta)*(number_of_recorded_dates-number_of_errors) ### not incorporating the binomial coefficient as we knoe exactly which ones are with and without error
  
  return(result)
}
# system.time(LL_error_term(aug_dat, theta, obs_dat))

#LL_error_term_slow<-function(aug_dat, theta, obs_dat)
#{
#  LL <- sum (sapply(1:n_groups, function(g) sum (sapply(1:ncol(aug_dat$D[[g]]), function(j) sum(LL_error_term_by_group_delay_and_indiv(aug_dat, theta, obs_dat, g, j, 1:nrow(obs_dat[[g]]))) ) ) ) )
#  return(LL)
#}
# system.time(LL_error_term_slow(aug_dat, theta, obs_dat))

#' @import EpiEstim
DiscrSI_vectorised <- function(x, mu, sigma, log=TRUE)
{
  if(log) res <- sapply(x, function(k) log(DiscrSI(k, mu+1, sigma))) else res <- sapply(x, function(k) DiscrSI(k, mu+1, sigma)) ### here we use mu+1 because we don't want the shifted gamma, just the gamma
  return(res)
}

DiscrSI_vectorised_from_mu_CV <- function(x, mu, CV, log=TRUE)
{
  sigma <- mu*CV
  if(log) res <- sapply(x, function(k) log(DiscrSI(k, mu+1, sigma))) else res <- sapply(x, function(k) DiscrSI(k, mu+1, sigma)) ### here we use mu+1 because we don't want the shifted gamma, just the gamma
  return(res)
}

LL_delays_term_by_group_delay_and_indiv <- function(aug_dat, theta, obs_dat, group_idx, delay_idx, indiv_idx, Delta=NULL)
{
  if(is.null(Delta)) Delta <- compute_delta_group_delay_and_indiv(aug_dat$D, group_idx, delay_idx, indiv_idx, index = index_dates)
  LL <- DiscrSI_vectorised_from_mu_CV(Delta + 1, theta$mu[[group_idx]][delay_idx], theta$CV[[group_idx]][delay_idx], log=TRUE)
  return(LL)
}

LL_delays_term<-function(aug_dat, theta, obs_dat, Delta=NULL)
{
  if(is.null(Delta)) Delta <- compute_delta(aug_dat$D, index = index_dates)
  LL <- sum (sapply(1:length(obs_dat), function(g) sum (sapply(2:ncol(aug_dat$D[[g]]), function(j) sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, theta, obs_dat, g, j-1, 1:nrow(obs_dat[[g]]), Delta[[g]][1:nrow(obs_dat[[g]]), j-1])) ) ) ) )
  return(LL)
}
# LL_delays_term(aug_dat, theta, obs_dat)

LL_total <- function(aug_dat, theta, obs_dat, range_dates=NULL)
{
  res <- LL_observation_term(aug_dat, theta, obs_dat, range_dates) + 
    LL_error_term(aug_dat, theta, obs_dat) + 
    LL_delays_term(aug_dat, theta, obs_dat)
  return(res)
}
# LL_total(aug_dat, theta, obs_dat)

###############################################
### priors ###
###############################################

# zeta ~ beta(low mean) # need a very informative prior
lprior_prob_error <- function(theta, shape1=3, shape2=12) 
{
  # can use this code to plot the corresponding prior: 
  # x <- seq(0,1,0.01)
  # y <- dbeta(x, shape1, shape2)
  # plot(x, y, type="l")
  return(dbeta(theta$zeta, shape1, shape2, log = TRUE))
}
# param_beta <- find_params_beta(mean=0.2, var=0.01)
# lprior_prob_error(theta, param_beta[1], param_beta[2])

# mu and CV ~ Exp(mean 1000) # very informative prior should be ok because data will be informative
lprior_params_delay <- function(what=c("mu", "CV"), theta, mean=100) # using the same prior for the mean of all delays
{
  what <- match.arg(what)
  # can use this code to plot the corresponding prior: 
  # x <- seq(0,1000,1)
  # y <- dexp(x, 1/mean)
  # plot(x, y, type="l")
  return(sum(dexp(unlist(theta[[what]]), 1/mean, log = TRUE)))
}
#lprior_params_delay("mu", theta)

lprior_total <- function(theta, shape1_prob_error=3, shape2_prob_error=12, mean_mean_delay=100, mean_std_delay=100)
{
  res <- lprior_prob_error(theta, shape1_prob_error, shape2_prob_error) + 
    lprior_params_delay("mu", theta, mean_mean_delay) + 
    lprior_params_delay("CV", theta, mean_mean_delay)
  return(res)
}
#lprior_total(theta)

###############################################
### posteriors ###
###############################################

lposterior_total <- function(aug_dat, theta, obs_dat, prior_shape1_prob_error=3, prior_shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100, range_dates=NULL)
{
  res <- LL_total(aug_dat, theta, obs_dat, range_dates) + lprior_total(theta, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay)
  return(res)
}
#lposterior_total(aug_dat, theta, obs_dat)