###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

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
  mu[[g]] <- rep(1.0,n_dates[g]-1)
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
        min_non_NA_value <- min(which(!is.na(D[[g]][e,])))-1
        for(f in min_non_NA_value:1)
        {
          D[[g]][e,f] <- D[[g]][e,f+1]
        }
      }
      if(any(is.na(D[[g]][e,]))) # dealing with remaining missing values if any
      {
        tmp <- which(is.na(D[[g]][e,]))
        for(f in tmp)
        {
          D[[g]][e,f] <- D[[g]][e,f-1]
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
### likelihood function ###
###############################################

LL_observation_term<-function(aug_dat, theta, obs_dat, log=TRUE)
{
  return(1)
}

LL_error_term<-function(aug_dat, theta, obs_dat, log=TRUE)
{
  number_of_errors<-0
  number_of_recorded_dates<-0
  for(g in 1:n_groups)
  {
    number_of_errors<-number_of_errors+sum(aug_dat$E[[g]]==1)
    number_of_recorded_dates<-number_of_recorded_dates+sum(aug_dat$E[[g]] != -1)
  }
  result<-dbinom(number_of_errors,number_of_recorded_dates,theta$zeta,log=log)

  return(result)
}
# LL_error_term(aug_dat, theta, obs_dat, log=TRUE)

LL_delays_term<-function()
{
  return(1)
}

