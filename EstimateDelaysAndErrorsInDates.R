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
### data ###
###############################################

raw_dat<-readRDS("Dat.rds")

# head(raw_dat)

N <- nrow(raw_dat) # Number of cases

colDates <- grep("Date", names(raw_dat))
tmp <- split(raw_dat[colDates],raw_dat$Path)
# splitting dataset according to Path and removing NA date columns in each of these
# - should only remain dates that are relevant for each group
obs_dat <- lapply(tmp, function(x) x[,colSums(is.na(x))!=nrow(x)] )

n_dates <- sapply(obs_dat, ncol )

n_groups <- length(n_dates)

###############################################
### define parameters to be used for initialisation of the chain ###
###############################################

### mean and std of distribution of various delays, by group
mu <- list()
for(g in 1:n_groups) 
{
  mu[[g]] <- rep(15.0,n_dates[g]-1)
}
names(mu) <- names(n_dates)
sigma <- mu

### list of all parameters
theta <- list(zeta = 0.05, # zeta is the probability for a date to be misrecorded, conditional on being recorded (<-> Ei != - 1)
              # TODO:
              # could consider having zeta being type of date specific (e.g. more error on onset than death dates),
              # time specific and/or space specific
              mu = mu, # mean of gamma distributions used to characterise the various delays in different groups: mu[[g]][k] is the mean k^th delay in group g
              sigma = sigma) # sigma of gamma distributions used to characterise the various delays in different groups: sigma[[g]][k] is the std k^th delay in group g


###############################################
### define augmented data to be used for initialisation of the chain ###
###############################################

### D contains the unobserved true dates ###

D <- list() ######### RICH MAYBE IMPROVE THIS BIT OF SHIT CODE ########
for(g in 1:n_groups) 
{
  D[[g]] <- obs_dat[[g]]
  for(e in 1:nrow(D[[g]]))
  {
    if(any(is.na(D[[g]][e,])))
    {
      tmp <- which(is.na(D[[g]][e,]))
      if(1 %in% tmp) # dealing with missing values ahead of the series of dates
      {
        min_non_NA_value <- min(which(!is.na(D[[g]][e,])))
        for(f in (min_non_NA_value-1):1)
        {
          D[[g]][e,index_dates[[g]][,match(f, index_dates[[g]][1,])][1]] <- D[[g]][e,index_dates[[g]][,match(f, index_dates[[g]][1,])][2]]
        }
      }
      if(any(is.na(D[[g]][e,]))) # dealing with remaining missing values if any
      {
        tmp <- which(is.na(D[[g]][e,]))
        for(f in tmp)
        {
          D[[g]][e,index_dates[[g]][,match(f, index_dates[[g]][2,])][2]] <- D[[g]][e,index_dates[[g]][,match(f, index_dates[[g]][2,])][1]]
        }
      }
    }
  }
}
names(D) <- names(obs_dat)

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
    E[[g]][,j] <- rbinom(nrow(obs_dat[[g]]), 1, theta$zeta)
    E[[g]][is.na(obs_dat[[g]][,j]),j] <- -1
  }
  names(E[[g]]) <- names(obs_dat[[g]])
}
names(E) <- names(obs_dat)

### now update D to be different to obs_dat if E = 1

for(g in 1:n_groups) 
{
  with_error <- which(E[[g]]==1, arr.ind =TRUE)
  for(ii in 1: nrow(with_error))
  {
    D[[g]][with_error[ii,1], with_error[ii,2]] <- obs_dat[[g]][with_error[ii,1], with_error[ii,2]] + sample(c(-1,1), 1)
  }
}

aug_dat <- list(D = D,
                E = E)


###############################################
### compute_delta function to compute relevant delays based on index, which tells you which dates should be used for delayl calculation ###
###############################################

compute_delta <- function(aug_dat, index = index_dates )
{
  Delta <- list()
  for(g in 1:n_groups)
  {
    Delta[[g]] <- matrix(NA, nrow(aug_dat$D[[g]]), ncol(aug_dat$D[[g]])-1)
    for(j in 2: ncol(aug_dat$D[[g]]))
    {
      Delta[[g]][,j-1] <- aug_dat$D[[g]][,index[[g]][,j-1][2]] - aug_dat$D[[g]][,index[[g]][,j-1][1]]
    }
  }
  return(Delta)
}

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

LL_observation_term<-function(aug_dat, theta, obs_dat)
{
  LL_no_error <- 0
  LL_error_or_missing <- 0
  range_dates <- find_range(obs_dat)
  for(g in 1:n_groups)
  {
    ### making sure D=y if E=0 ### note could remove this if by construction this is always true - could speed up code
    no_error <- aug_dat$E[[g]]==0
    LL_no_error <- LL_no_error + sum( log(aug_dat$D[[g]][no_error] == obs_dat[[g]][no_error]) )
    ### if E=1, what is the relationship between true date D and observed date y
    # for now, observation likelihood conditional on E=1 is uniform on the range of observed dates
    ### same for E=-1, D can take any value in range with same probability as they are all consistent with y=NA
    error_or_missing <- (aug_dat$E[[g]]==-1 | aug_dat$E[[g]]==1)
    ### K is the relative probability of observing a given error, conditional on presence of error
    # for now, K is given as 1/n, where n is the number of dates in the range_dates
    # could use something different if we define the space of possible errors differently. 
    K <- (1/as.numeric(diff(range_dates))) ### think about impact of this choice
    LL_error_or_missing <- LL_error_or_missing + sum( log( K* ((aug_dat$D[[g]][error_or_missing] >= range_dates[1]) & (aug_dat$D[[g]][error_or_missing] <= range_dates[2])) ) )
  }
  LL <- LL_no_error + LL_error_or_missing
  if(is.infinite(LL)) LL <- -100000 # arbitrarily small number to avoid -Inf
  
  return(LL)
}
# LL_observation_term(aug_dat, theta, obs_dat)

LL_error_term<-function(aug_dat, theta, obs_dat)
{
  number_of_errors<-0
  number_of_recorded_dates<-0
  for(g in 1:n_groups)
  {
    number_of_errors<-number_of_errors+sum(aug_dat$E[[g]]==1)
    number_of_recorded_dates<-number_of_recorded_dates+sum(aug_dat$E[[g]] != -1)
  }
  result<-dbinom(number_of_errors,number_of_recorded_dates,theta$zeta,log=TRUE)
  
  return(result)
}
# LL_error_term(aug_dat, theta, obs_dat)

DiscrSI_vectorised <- function(x, mu, sigma, log=TRUE)
{
  if(log) res <- sapply(x, function(k) log(DiscrSI(k, mu, sigma))) else res <- sapply(x, function(k) DiscrSI(k, mu, sigma))
  return(res)
}

LL_delays_term<-function(aug_dat, theta, obs_dat)
{
  Delta <- compute_delta(aug_dat)
  LL <- 0
  for(g in 1:n_groups)
  {
    for(j in 2:ncol(aug_dat$D[[g]]))
    {
      LL_vect_this_group <- DiscrSI_vectorised(Delta[[g]] + 1, theta$mu[[g]][j-1], theta$sigma[[g]][j-1], log=TRUE)
      LL <- LL + sum( LL_vect_this_group )
    }
  }
  return(LL)
}
# LL_delays_term(aug_dat, theta, obs_dat)

LL_total <- function(aug_dat, theta, obs_dat)
{
  res <- LL_observation_term(aug_dat, theta, obs_dat) + 
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
lprior_prob_error <- function(theta, mean=0.2, var=0.01) 
{
  param_beta <- find_params_beta(mean, var)
  # can use this code to plot the corresponding prior: 
  # x <- seq(0,1,0.01)
  # y <- dbeta(x, param_beta[1], param_beta[2])
  # plot(x, y, type="l")
  return(dbeta(theta$zeta, param_beta[1], param_beta[2], log = TRUE))
}
#lprior_prob_error(theta)

# mu and CV ~ Exp(mean 1000) # very informative prior should be ok because data will be informative
lprior_mean_delay <- function(theta, mean=100) # using the same prior for the mean of all delays
{
  # can use this code to plot the corresponding prior: 
  # x <- seq(0,1000,1)
  # y <- dexp(x, 1/mean)
  # plot(x, y, type="l")
  return(sum(dexp(unlist(theta$mu), 1/mean, log = TRUE)))
}
#lprior_mean_delay(theta)

### NEED TO CHANGE THIS TO BE THE CV RATHER THAN STD
lprior_std_delay <- function(theta, mean=100) # using the same prior for the std of all delays
{
  # can use this code to plot the corresponding prior: 
  # x <- seq(0,1000,1)
  # y <- dexp(x, 1/mean)
  # plot(x, y, type="l")
  return(sum(dexp(unlist(theta$sigma), 1/mean, log = TRUE)))
}
#lprior_std_delay(theta)

lprior_total <- function(theta, mean_prob_error=0.2, var_prob_error=0.01, mean_mean_delay=100, mean_std_delay=100)
{
  res <- lprior_prob_error(theta, mean_prob_error, var_prob_error) + 
    lprior_mean_delay(theta, mean_mean_delay) + 
    lprior_std_delay(theta, mean_mean_delay)
  return(res)
}
#lprior_total(theta)

###############################################
### posteriors ###
###############################################

### TO WRITE

###############################################
### move functions ###
###############################################

### TO WRITE

### plan: 

## D_i ##
# probability of moving or not
# moves by +/-1 random walk, and E_i is adjusted accordingly , 
# i.e. if D_i moves to y_i then E_i moves to 0, else E_i moves to 1. 
# then accept/reject based on posterior values
# this is symmetrical so no correction needed

## alpha, beta ## --> propose to reparameterise to mean and CV
# move with a lognormal proposal, with adequate correction

## zeta ##
# move with a truncated (<1) lognormal proposal, with adequate correction

###############################################
### MCMC ###
###############################################

### TO WRITE

###############################################
### TO DO ###
###############################################

# Anne: 
# find a good starting point for aug_data
# fill in bits of codes already existing to various places

# Marc: 
# finish writing
# think about the 1/(T-T0)

# Future ideas: 
# simulation study
# other datasets - Marc to talk to John? 
# outputs: proportion erroneous data - ... - nice graphs




