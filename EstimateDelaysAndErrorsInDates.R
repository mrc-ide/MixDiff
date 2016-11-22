###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

library(EpiEstim) # to use DiscrSI which does the discretised Gamma

###############################################
### index_dates says which combinations of dates to use for the delays ###
###############################################

index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )

###############################################
### index_dates_order says which dates should be before which ones - used for initialisation ###
###############################################

index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )

###############################################
### function to handle dates ###
###############################################

date_to_int <- function(date, origin = "1970-01-01")
{
  return(as.integer(date - as.Date(origin)))
}

int_to_date <- function(int, origin = "1970-01-01")
{
  return(int + as.Date(origin))
}

###############################################
### data ###
###############################################

raw_dat<-readRDS("Dat.rds")

# head(raw_dat)

N <- nrow(raw_dat) # Number of cases

colDates <- grep("Date", names(raw_dat))
tmp <- split(raw_dat[colDates],raw_dat$Path)
# splitting dataset according to Path and removing NA date columns in each of these
# - should only remain dates that are relevant for each group
obs_dat <- lapply(tmp, function(x) sapply(which(colSums(is.na(x))!=nrow(x)), function(j) date_to_int(x[,j]) )) ### converting obs_dat to be integers - easier to handle than dates

n_dates <- sapply(obs_dat, ncol )

n_groups <- length(n_dates)

###############################################
### compute_delta function to compute relevant delays based on index, which tells you which dates should be used for delayl calculation ###
###############################################

compute_delta_group_delay_and_indiv<-function(D, group_idx, delay_idx, indiv_idx, index = index_dates)
{
  Delta <- D[[group_idx]][indiv_idx,index[[group_idx]][,delay_idx][2]] - D[[group_idx]][indiv_idx,index[[group_idx]][,delay_idx][1]]
  return(Delta)
}

compute_delta <- function(D, index = index_dates)
{
  Delta <- lapply(1:n_groups, function(g){
    m <- matrix(NA, nrow(D[[g]]), ncol(D[[g]])-1)
    for(j in 1: ncol(m))
    {
      m[,j] <- D[[g]][,index[[g]][,j][2]] - D[[g]][,index[[g]][,j][1]]
    }
    return(m)
  })
  return(Delta)
}

###############################################
### define parameters to be used for initialisation of the chain ###
###############################################

### mean and std of distribution of various delays, by group
### we use a the starting point the observed mean and std of each delay in each group
obs_delta <- compute_delta(obs_dat, index_dates)
mu <- lapply(1:n_groups, function(g) abs(apply(obs_delta[[g]], 2, mean, na.rm=TRUE) ))
names(mu) <- names(n_dates)
sigma <- lapply(1:n_groups, function(g) abs(apply(obs_delta[[g]], 2, sd, na.rm=TRUE) ))
names(sigma) <- names(n_dates)

### list of all parameters
theta <- list(zeta = 0.01, # zeta is the probability for a date to be misrecorded, conditional on being recorded (<-> Ei != - 1)
              # TODO:
              # could consider having zeta being type of date specific (e.g. more error on onset than death dates),
              # time specific and/or space specific
              mu = mu, # mean of gamma distributions used to characterise the various delays in different groups: mu[[g]][k] is the mean k^th delay in group g
              sigma = sigma) # sigma of gamma distributions used to characterise the various delays in different groups: sigma[[g]][k] is the std k^th delay in group g


###############################################
### define augmented data to be used for initialisation of the chain ###
###############################################

### D contains the unobserved true dates ###

initialise_aug_data <- function(obs_dat, index_dates_order)
{
  D <- list()
  for(g in 1:n_groups) 
  {
    D[[g]] <- obs_dat[[g]]
    for(e in 1:nrow(D[[g]]))
    {
      #print(e)
      missing_dates <- which(is.na(D[[g]][e,]))
      while(length(missing_dates)>0)
      {
        can_be_inferred_from <- lapply(missing_dates, function(i) {
          x <- which(index_dates_order[[g]]==i, arr.ind = TRUE)
          from_idx <- sapply(1:nrow(x), function(k) index_dates_order[[g]][-x[k,1],x[k,2]] )
          from_value <- sapply(1:nrow(x), function(k) D[[g]][e,index_dates_order[[g]][-x[k,1],x[k,2]]])
          rule <- sapply(1:nrow(x), function(k) if(x[k,1]==1) "before" else "after"  )
          return(list(rule=rule,from_idx=from_idx, from_value=from_value))
        })
        can_be_inferred <- which(sapply(1:length(missing_dates), function(i) any(!is.na(can_be_inferred_from[[i]]$from_value))))
        for(k in can_be_inferred)
        {
          x <- which(!is.na(can_be_inferred_from[[k]]$from_value))
          if(length(x)==1)
          {
            inferred <- can_be_inferred_from[[k]]$from_value[x]
          }else
          {
            if(all(can_be_inferred_from[[k]]$rule[x] == "before"))
            {
              inferred <- min(can_be_inferred_from[[k]]$from_value[x])
            }else if (all(can_be_inferred_from[[k]]$rule[x] == "after"))
            {
              inferred <- max(can_be_inferred_from[[k]]$from_value[x])
            }else
            {
              max_val <- min(can_be_inferred_from[[k]]$from_value[x][can_be_inferred_from[[k]]$rule[x] %in% "before"])
              min_val <- max(can_be_inferred_from[[k]]$from_value[x][can_be_inferred_from[[k]]$rule[x] %in% "after"])
              if(min_val>max_val)
              {
                stop("Incompatible data to infer from. ")
              }else
              {
                inferred <- median(c(min_val, max_val))
              }
            }
          }
          D[[g]][e,missing_dates[k]] <- inferred
        }
        missing_dates <- which(is.na(D[[g]][e,]))
      }
    }
  }
  names(D) <- names(obs_dat)
  return(D)
}

D <- initialise_aug_data(obs_dat, index_dates_order)

### E contains an indicator of whether the observed date is the true one or not: ###
### E = -1 if date is unobserved i.e. obs_dat = -1 ###
### E = 0 if date is observed and exact i.e. obs_dat = D ###
### E = 1 if date is observed and unexact i.e. obs_dat not necessarily = D ###

E <- list()
for(g in 1:n_groups) 
{
  E[[g]] <- as.data.frame(matrix(NA,nrow(obs_dat[[g]]),ncol(obs_dat[[g]])))
  for(j in 1:ncol(obs_dat[[g]]))
  {
    E[[g]][,j] <- 0 # assume no error to start with ### rbinom(nrow(obs_dat[[g]]), 1, theta$zeta)
    E[[g]][is.na(obs_dat[[g]][,j]),j] <- -1
  }
  names(E[[g]]) <- names(obs_dat[[g]])
}
names(E) <- names(obs_dat)

### now update D to be different to obs_dat if E = 1
# for(g in 1:n_groups) 
# {
#   with_error <- which(E[[g]]==1, arr.ind =TRUE)
#   for(ii in 1: nrow(with_error))
#   {
#     D[[g]][with_error[ii,1], with_error[ii,2]] <- obs_dat[[g]][with_error[ii,1], with_error[ii,2]] + sample(c(-1,1), 1)
#   }
# }

aug_dat <- list(D = D,
                E = E)

###############################################
### likelihood function ###
###############################################

find_range <- function(obs_dat)
{
  min_date <- min(obs_dat[[1]][,1], na.rm=TRUE)
  max_date <- max(obs_dat[[1]][,1], na.rm=TRUE)
  for(g in 1:n_groups)
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
  LL <- sum (sapply(1:n_groups, function(g) sum (sapply(1:ncol(aug_dat$D[[g]]), function(j) sum(LL_observation_term_by_group_delay_and_indiv(aug_dat, theta, obs_dat, g, j, 1:nrow(obs_dat[[g]]),range_dates)) ) ) ) )
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

DiscrSI_vectorised <- function(x, mu, sigma, log=TRUE)
{
  if(log) res <- sapply(x, function(k) log(DiscrSI(k, mu, sigma))) else res <- sapply(x, function(k) DiscrSI(k, mu, sigma))
  return(res)
}

LL_delays_term_by_group_delay_and_indiv <- function(aug_dat, theta, obs_dat, group_idx, delay_idx, indiv_idx, Delta=NULL)
{
  if(is.null(Delta)) Delta <- compute_delta_group_delay_and_indiv(aug_dat$D, group_idx, delay_idx, indiv_idx, index = index_dates)
  LL <- DiscrSI_vectorised(Delta + 1, theta$mu[[group_idx]][delay_idx], theta$sigma[[group_idx]][delay_idx], log=TRUE)
  return(LL)
}

LL_delays_term<-function(aug_dat, theta, obs_dat, Delta=NULL)
{
  if(is.null(Delta)) Delta <- compute_delta(aug_dat$D, index = index_dates)
  LL <- sum (sapply(1:n_groups, function(g) sum (sapply(2:ncol(aug_dat$D[[g]]), function(j) sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, theta, obs_dat, g, j-1, 1:nrow(obs_dat[[g]]), Delta[[g]][1:nrow(obs_dat[[g]]), j-1])) ) ) ) )
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
# test if works: 
# param_beta <- find_params_beta(0.1, 0.05)
# sample_beta <- rbeta(1000, param_beta[1], param_beta[2])
# mean(sample_beta)
# var(sample_beta)

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
lprior_params_delay <- function(what=c("mu", "sigma"), theta, mean=100) # using the same prior for the mean of all delays
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
    lprior_params_delay("sigma", theta, mean_mean_delay)
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

###############################################
### move functions ###
###############################################

## D_i ##
# probability of moving or not
# moves by +/-1 random walk, and E_i is adjusted accordingly , 
# i.e. if D_i moves to y_i then E_i moves to 0, else E_i moves to 1. 
# then accept/reject based on posterior values
# this is symmetrical so no correction needed
move_Di <- function(i, group_idx, date_idx, 
                    curr_aug_dat,
                    theta, 
                    obs_dat, 
                    shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100, range_dates=NULL) 
{
  if(is.null(range_dates)) range_dates <- find_range(obs_dat)
  
  # draw proposed value for D using +/-1 random walk
  curr_aug_dat_value <- curr_aug_dat$D[[group_idx]][i,date_idx]
  proposed_aug_dat_value <- curr_aug_dat_value + sample(c(-1,1), length(i), replace=TRUE)
  
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
  # ratio_post_long <- lposterior_total(proposed_aug_dat, theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay) - 
  # lposterior_total(curr_aug_dat, theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay)
  
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
# test_move_Di <- move_Di(i=1, group_idx=1, date_idx=1, curr_aug_dat = aug_dat, theta, obs_dat, shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100) 
# test_move_Di$new_aug_dat$D[[1]][1,1] # new value
# aug_dat$D[[1]][1,1] # old value

### move mu with a lognormal proposal ###   # NOTE: consider changing sigma to be CV
move_lognormal <- function(what=c("mu","sigma"), group_idx, delay_idx, sdlog, 
                           aug_dat,
                           curr_theta, 
                           obs_dat, 
                           shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100) 
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
  }else if(what=="sigma")
  {
    ratio_post <- lprior_params_delay(what, proposed_theta, prior_mean_std_delay) - lprior_params_delay(what, curr_theta, prior_mean_std_delay) 
  }
  Delta <- compute_delta_group_delay_and_indiv(aug_dat$D, group_idx, delay_idx,  1:nrow(obs_dat[[group_idx]]), index = index_dates) # same for proposed and curent par values so no need to recompute twice
  ratio_post <- ratio_post + sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, proposed_theta, obs_dat, group_idx, delay_idx, 1:nrow(obs_dat[[group_idx]])), Delta) - 
    sum(LL_delays_term_by_group_delay_and_indiv(aug_dat, curr_theta, obs_dat, group_idx, delay_idx, 1:nrow(obs_dat[[group_idx]])), Delta) 
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(aug_dat, proposed_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay) - 
  # lposterior_total(aug_dat, curr_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay)
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
# test_move_mu <- move_lognormal(what="mu", group_idx=1, delay_idx=1, sdlog=0.1, aug_dat, curr_theta = theta, obs_dat, shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100)
# test_move_mu$new_theta$mu[[1]][1] # new value
# theta$mu[[1]][1] # old value

### move zeta with a truncated (<1) lognormal proposal ### NOT USED ANYMORE, NOW USING GIBBS SAMPLE FOR ZETA
move_truncated_lognormal <- function(what=c("zeta"), sdlog, upper_bound=1,  
                                     aug_dat,
                                     curr_theta, 
                                     obs_dat, 
                                     shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100) 
{
  what <- match.arg(what)
  
  # draw proposed value
  curr_param_value <- curr_theta[[what]]
  proposed_param_value <- rlnorm(1,meanlog=log(curr_param_value), sdlog=sdlog)
  while(proposed_param_value > upper_bound) # repeat until draw a value below the upper bound
  {
    proposed_param_value <- rlnorm(1,meanlog=log(curr_param_value), sdlog=sdlog)
  }
  
  proposed_theta <- curr_theta
  proposed_theta[[what]] <- proposed_param_value
  
  # calculates probability of acceptance
  ratio_post <- lprior_prob_error(proposed_theta, shape1=shape1_prob_error, shape2=shape2_prob_error) - 
    lprior_prob_error(curr_theta, shape1=shape1_prob_error, shape2=shape2_prob_error)
  ratio_post <- ratio_post + LL_error_term(aug_dat, proposed_theta, obs_dat) -
    LL_error_term(aug_dat, curr_theta, obs_dat) 
  
  ### note that ratio_post should be the same as: 
  # ratio_post_long <- lposterior_total(aug_dat, proposed_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay) - 
  # lposterior_total(aug_dat, curr_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay)
  
  correction <- log(proposed_param_value) - log(curr_param_value) # correction for lognormal distribution
  tmp1 <- ( log(upper_bound) - log(curr_param_value) ) / sdlog
  tmp2 <- ( log(upper_bound) - log(proposed_param_value) ) / sdlog
  correction <- correction + pnorm(tmp1, log.p = TRUE) - pnorm(tmp2, log.p = TRUE)  # additional correction for the truncation
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
# test_move_zeta <- move_truncated_lognormal(what="zeta", sdlog=0.005, upper_bound=1,  aug_dat, curr_theta = theta, obs_dat, shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100)
# test_move_zeta$new_theta$zeta # new value
# theta$zeta # old value

### move zeta with a truncated (<1) lognormal proposal ###
move_zeta_gibbs <- function(aug_dat,
                            curr_theta, 
                            obs_dat, 
                            shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100) 
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
# test_move_zeta_gibbs <- move_zeta_gibbs(aug_dat, curr_theta = theta, obs_dat, shape1_prob_error=3, shape2_prob_error=12, prior_mean_mean_delay=100, prior_mean_std_delay=100)
# test_move_zeta_gibbs$new_theta$zeta # new value
# theta$zeta # old value

###############################################
### MCMC ###
###############################################

n_iter <- 10 # currently (21st Nov 2016, updating 1/10th of Di per group at each iteration, 100 iterations take ~360 seconds)

move_D_by_groups_of_size <- 1

### prior parameters 

prior_shape1_prob_error=3
prior_shape2_prob_error=12
prior_mean_mean_delay=100
prior_mean_std_delay=100

### initialisation

range_dates <- find_range(obs_dat)

# to store param values
curr_theta <- theta
theta_chain <- curr_theta

# to store augmented data values
curr_aug_dat <- aug_dat
aug_dat_chain <- list(D=list(), E=list())
for(g in 1:n_groups)
{
  aug_dat_chain$D[[g]] <- lapply(1:n_dates[[g]], function(j) as.integer(aug_dat$D[[g]][,j]))
  names(aug_dat_chain$D[[g]]) <- paste0("Delay",1:n_dates[[g]])
  aug_dat_chain$E[[g]] <- lapply(1:n_dates[[g]], function(j) as.integer(aug_dat$E[[g]][,j]))
  names(aug_dat_chain$E[[g]]) <- paste0("Delay",1:n_dates[[g]])
}
names(aug_dat_chain$D) <- names(obs_dat)
names(aug_dat_chain$E) <- names(obs_dat)

add_new_value_chain_theta <- function(theta_chain, new_theta)
{
  theta_chain$zeta <- c(theta_chain$zeta, new_theta$zeta)
  for(g in 1:n_groups)
  {
    theta_chain$mu[[g]] <- rbind(theta_chain$mu[[g]], new_theta$mu[[g]])
    theta_chain$sigma[[g]] <- rbind(theta_chain$sigma[[g]], new_theta$sigma[[g]])
  }
  return(theta_chain)
}

add_new_value_chain_aug_dat <- function(curr_aug_dat, new_aug_dat)
{
  for(g in 1:n_groups)
  {
    for(j in 1:n_dates[[g]])
    {
      curr_aug_dat$D[[g]][[j]] <- rbind(curr_aug_dat$D[[g]][[j]], new_aug_dat$D[[g]][,j])
      curr_aug_dat$E[[g]][[j]] <- rbind(curr_aug_dat$E[[g]][[j]], new_aug_dat$E[[g]][,j])
    }
  }
  return(curr_aug_dat)
}

logpost_chain <- rep(NA, n_iter)
logpost_chain[1] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay, range_dates)

n_accepted_D_moves <- 0
n_proposed_D_moves <- 0

n_accepted_mu_moves <- 0
n_proposed_mu_moves <- 0

n_accepted_sigma_moves <- 0
n_proposed_sigma_moves <- 0

#n_accepted_zeta_moves <- 0 # not used as Gibbs sampler 
#n_proposed_zeta_moves <- 0 # not used as Gibbs sampler

### turn on and off various moves
D_moves_on <- TRUE
mu_moves_on <- TRUE
sigma_moves_on <- TRUE
zeta_moves_on <- TRUE

### std of moves

fraction_Di_to_update <- 1/10

sdlog_mu <- 0.1 # for now moving all mus with the same sd, 
# might need to revisit this as some delays might be longer than others an require different sdlog to optimise mixing of the chain

sdlog_sigma <- 0.05 # for now moving all sigmas with the same sd, 
# might need to revisit this as some delays might be longer than others an require different sdlog to optimise mixing of the chain

#sdlog_zeta <- 0.005 # not used as Gibbs sampler

system.time({
  for(k in 1:(n_iter-1))
  {
    print(k)
    
    # move some of the D_i
    if(D_moves_on)
    {
      for(g in 1:n_groups)
      {
        for(j in 1:ncol(curr_aug_dat$D[[g]]))
        {
          to_update <- sample(1:nrow(obs_dat[[g]]), round(nrow(obs_dat[[g]])*fraction_Di_to_update)) # proposing moves for only a certain fraction of dates
          n_10_to_update <- floor(length(to_update) / move_D_by_groups_of_size)
          for(i in 1:length(n_10_to_update))
          {
            tmp <- move_Di (to_update[move_D_by_groups_of_size*(i-1)+(1:move_D_by_groups_of_size)], g, j, 
                            curr_aug_dat,
                            curr_theta, 
                            obs_dat, 
                            prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay, range_dates) 
            n_proposed_D_moves <- n_proposed_D_moves + 1
            n_accepted_D_moves <- n_accepted_D_moves + tmp$accept
            if(tmp$accept==1) curr_aug_dat <- tmp$new_aug_dat # if accepted move, update accordingly
          }
        }
      }
    }
    
    # move zeta using Gibbs sampler
    if(zeta_moves_on)
    {
      tmp <- move_zeta_gibbs(curr_aug_dat,
                             curr_theta, 
                             obs_dat, 
                             prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay) 
      curr_theta <- tmp$new_theta # always update with new theta (Gibbs sampler)
    }
    
    # move mu
    if(mu_moves_on)
    {
      for(g in 1:n_groups)
      {
        for(j in 2:ncol(curr_aug_dat$D[[g]]))
        {
          tmp <- move_lognormal(what="mu", g, j-1, sdlog_mu, 
                                curr_aug_dat,
                                curr_theta, 
                                obs_dat, 
                                prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay)
          n_proposed_mu_moves <- n_proposed_mu_moves + 1
          n_accepted_mu_moves <- n_accepted_mu_moves + tmp$accept
          if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
        }
      }
    }
    
    # move sigma
    if(sigma_moves_on)
    {
      for(g in 1:n_groups)
      {
        for(j in 2:ncol(curr_aug_dat$D[[g]]))
        {
          tmp <- move_lognormal(what="sigma", g, j-1, sdlog_sigma, 
                                curr_aug_dat,
                                curr_theta, 
                                obs_dat, 
                                prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay)
          n_proposed_sigma_moves <- n_proposed_sigma_moves + 1
          n_accepted_sigma_moves <- n_accepted_sigma_moves + tmp$accept
          if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
        }
      }
    }
    
    # recording the value of parameters after all moves
    theta_chain <- add_new_value_chain_theta(theta_chain, curr_theta)
    aug_dat_chain <- add_new_value_chain_aug_dat(aug_dat_chain, curr_aug_dat)
    
    # recording the likelihood after all moves
    logpost_chain[k+1] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_std_delay, range_dates)
  }
})

# save.image("tmp.Rdata")

###############################################
### acceptance probabilities ###
###############################################

n_accepted_D_moves / n_proposed_D_moves
#n_accepted_zeta_moves / n_proposed_zeta_moves # not computed as now using Gibbs sampler for Zeta
n_accepted_mu_moves / n_proposed_mu_moves
n_accepted_sigma_moves / n_proposed_sigma_moves

###############################################
### remove burnin ###
###############################################

burnin <- 1 # 1:100
logpost_chain <- logpost_chain[-burnin]
theta_chain$zeta <- theta_chain$zeta[-burnin]
for(g in 1:n_groups)
{
  theta_chain$mu[[g]] <- theta_chain$mu[[g]][-burnin,]
  theta_chain$sigma[[g]] <- theta_chain$sigma[[g]][-burnin,]
}
for(g in 1:n_groups)
{
  for(j in 1:n_dates[[g]])
  {
    aug_dat_chain$D[[g]][[j]] <- aug_dat_chain$D[[g]][[j]][-burnin,]
    aug_dat_chain$E[[g]][[j]] <- aug_dat_chain$E[[g]][[j]][-burnin,]
  }
}

###############################################
### plotting the MCMC output ###
###############################################

### parameters ###

pdf("ParamConvergencePlots.pdf", width=14, height=7)
par(mfrow=c(2, 5),mar=c(5, 6, 1, 1))

# looking at the logposterior chain 
plot(logpost_chain, type="l", xlab="Iterations", ylab="Log posterior")

# looking at mean delay 
group_idx <- 1 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(non hospitalised-alive group)", ylim=c(0, 15))
legend("topright", "Onset-Report", lty=1)
group_idx <- 2 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]][,j]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(non hospitalised-dead group)", ylim=c(0, 15))
for(j in 2:(n_dates[group_idx]-1))
{
  mu <- theta_chain$mu[[group_idx]][,j]
  lines(mu, col=j)
}
legend("topright", c("Onset-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 3 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]][,j]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(hospitalised-alive group)", ylim=c(0, 15))
for(j in 2:(n_dates[group_idx]-1))
{
  mu <- theta_chain$mu[[group_idx]][,j]
  lines(mu, col=j)
}
legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 4 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]][,j]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(hospitalised-dead group)", ylim=c(0, 15))
for(j in 2:(n_dates[group_idx]-1))
{
  mu <- theta_chain$mu[[group_idx]][,j]
  lines(mu, col=j)
}
legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])

# looking at zeta
zeta <- theta_chain$zeta
plot(zeta, type="l", xlab="Iterations", ylab="zeta")

# looking at std delay
group_idx <- 1 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(non hospitalised-alive group)", ylim=c(0, 15))
legend("topright", "Onset-Report", lty=1)
group_idx <- 2 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]][,j]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(non hospitalised-dead group)", ylim=c(0, 15))
for(j in 2:(n_dates[group_idx]-1))
{
  sigma <- theta_chain$sigma[[group_idx]][,j]
  lines(sigma, col=j)
}
legend("topright", c("Onset-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 3 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]][,j]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(hospitalised-alive group)", ylim=c(0, 15))
for(j in 2:(n_dates[group_idx]-1))
{
  sigma <- theta_chain$sigma[[group_idx]][,j]
  lines(sigma, col=j)
}
legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 4 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]][,j]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(hospitalised-dead group)", ylim=c(0, 15))
for(j in 2:(n_dates[group_idx]-1))
{
  sigma <- theta_chain$sigma[[group_idx]][,j]
  lines(sigma, col=j)
}
legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
dev.off()

### augmented data ###

pdf("AugDataConvergencePlots.pdf", width=14, height=14)
par(mfrow=c(4, 5),mar=c(5, 6, 1, 1))
group_idx <- 1 ##########################
# randomly pick 5 individuals in that group
indiv_to_plot <- sample(1:ncol(aug_dat_chain$D[[group_idx]][[1]]), 5)
for(i in 1:length(indiv_to_plot)) 
{
  j <- 1
  date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
  plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
  for(j in 2:(n_dates[group_idx]))
  {
    date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
    lines(date, col=j)
  }
  legend("topright", c("Onset","Report"), lty=1, col=1:n_dates[group_idx])
}
group_idx <- 2 ##########################
# randomly pick 5 individuals in that group
indiv_to_plot <- sample(1:ncol(aug_dat_chain$D[[group_idx]][[1]]), 5)
for(i in 1:length(indiv_to_plot)) 
{
  j <- 1
  date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
  plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
  for(j in 2:(n_dates[group_idx]))
  {
    date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
    lines(date, col=j)
  }
  legend("topright", c("Onset","Death","Report"), lty=1, col=1:n_dates[group_idx])
}
group_idx <- 3 ##########################
# randomly pick 5 individuals in that group
indiv_to_plot <- sample(1:ncol(aug_dat_chain$D[[group_idx]][[1]]), 5)
for(i in 1:length(indiv_to_plot)) 
{
  j <- 1
  date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
  plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
  for(j in 2:(n_dates[group_idx]))
  {
    date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
    lines(date, col=j)
  }
  legend("topright", c("Onset","Hosp","Disch","Report"), lty=1, col=1:n_dates[group_idx])
}
group_idx <- 4 ##########################
# randomly pick 5 individuals in that group
indiv_to_plot <- sample(1:ncol(aug_dat_chain$D[[group_idx]][[1]]), 5)
for(i in 1:length(indiv_to_plot)) 
{
  j <- 1
  date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
  plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
  for(j in 2:(n_dates[group_idx]))
  {
    date <- aug_dat_chain$D[[group_idx]][[j]][,indiv_to_plot[i]]
    lines(date, col=j)
  }
  legend("topright", c("Onset","Hosp","Death","Report"), lty=1, col=1:n_dates[group_idx])
}
dev.off()


### example of one that mixes better: 
#plot(aug_dat_chain$D[[1]][[2]][,6], type="l")

### examining the number of accepted changes per group and date
find_number_successful_changes <- function(aug_dat_chain, group_idx, date_idx)
{
  tmp <- abs(apply(aug_dat_chain$D[[group_idx]][[date_idx]], 2, diff))
  ret <- table(colSums(tmp)) # this says how many dates have had no changes, 1 change, etc... for that group and that date
  # cbind(colSums(tmp), obs_dat[[group_idx]]) ## the ones that are moved successfully are the ones with missing data
  return(ret)
}
number_successful_changes <- lapply(1:n_groups, function(group_idx) lapply(1:n_dates[group_idx], function(date_idx) find_number_successful_changes(aug_dat_chain, group_idx, date_idx)))
# so most of them seem stuck; they never move

sapply(1:n_groups, function(g) sapply(1:n_dates[g], function(j) unique(as.vector(aug_dat_chain$E[[g]][[j]])) ) )
# so we never reach E=1 i.e. there is an error. 
# This is because the probability of observing a given error is very very very small compared to the probability of observing no error
# with the range of observed dates, we have 
# K = 1/as.numeric(diff(range_dates)) # 0.001760563
# assuming zeta is 0.1, K*zeta = 0.0001760563 VS (1-zeta) = 0.9
# even worst for smaller values of zeta

###############################################
### TO DO ###
###############################################

# Anne: 
# check the MCMC, 
# try to speed up if possible
# considering only calculating the likelihood for some iterations (e.g. after burnin and thinning), posthoc? 
# should we update zeta after each D_i move, or after all D_i in a group move? 
# keep track of acceptance rate for D and for mu/sigma per group and per deay rather than altogether, to check if some moves are more successful than others. 
# also consider using Gibbs samplers to move mu and sigma --> for this need to reformulate as shape/scale: but doesn't seem obvious to sample from the posterior distribution? 

# Marc: 
# finish writing
# think about the 1/(T-T0)

# Future ideas: 
# simulation study
# other datasets - Marc to talk to John? 
# outputs: proportion erroneous data - ... - nice graphs




