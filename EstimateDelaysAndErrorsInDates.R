###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

###############################################
### data ###
###############################################

dat<-readRDS("Dat.rds")

# head(dat)

N <- nrow(dat) # Number of cases

colDates <- grep("Date", names(dat))
tmp <- split(dat[colDates],dat$Path)
# splitting dataset according to Path and removing NA date columns in each of these
# - should only remain dates that are relevant for each group
dat_by_group <- lapply(tmp, function(x) x[,colSums(is.na(x))!=nrow(x)] )

n_dates <- sapply(dat_by_group, ncol )

###############################################
### define parameters ###
###############################################

### mean and std of distribution of various delays, by group
mu <- list()
for(i in 1:length(n_dates))
{
  mu[[i]] <- rep(1.0,n_dates[i]-1)
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
### define augmented data ###
###############################################

### D contains the unobserved true dates ###

D <- list() ######### RICH MAYBE IMPROVE THIS BIT OF SHIT CODE ########
for(i in 1:length(n_dates))
{
  D[[i]] <- dat_by_group[[i]]
  for(e in 1:nrow(D[[i]]))
  {
    if(any(is.na(D[[i]][e,])))
    {
      tmp <- which(is.na(D[[i]][e,]))
      if(1 %in% tmp) # dealing with missing values ahead of the series of dates
      {
        min_non_NA_value <- min(which(!is.na(D[[i]][e,])))-1
        for(f in min_non_NA_value:1)
        {
          D[[i]][e,f] <- D[[i]][e,f+1]
        }
      }
      if(any(is.na(D[[i]][e,]))) # dealing with remaining missing values if any
      {
        tmp <- which(is.na(D[[i]][e,]))
        for(f in tmp)
        {
          D[[i]][e,f] <- D[[i]][e,f-1]
        }
      }
    }
  }
}
names(D) <- names(dat_by_group)

### E contains an indicator of whether the observed date is the true one or not: ###
### E = -1 if date is unobserved i.e. dat_by_group = -1 ###
### E = 0 if date is observed and exact i.e. dat_by_group = D ###
### E = 1 if date is observed and unexact i.e. dat_by_group not necessarily = D ###

E <- list()
for(i in 1:length(n_dates))
{
  E[[i]] <- as.data.frame(matrix(NA,nrow(dat_by_group[[i]]),ncol(dat_by_group[[i]])))
  for(j in 1:ncol(dat_by_group[[i]]))
  {
    E[[i]][,j] <- rbinom(nrow(dat_by_group[[i]]), 1, theta$zeta)
    E[[i]][is.na(dat_by_group[[i]][,j]),j] <- -1
  }
  names(E[[i]]) <- names(dat_by_group[[i]])
}
names(E) <- names(dat_by_group)

### Note that at the moment E can be 1 i.e. there is an error, yet D=dat_by_group i.e. in effect the error is null
### TO DO: TRY AND CHANGE SO THAT E=0 <-> dat_by_group=D

aug_dat <- list(D = D,
                E = E)


###############################################
### likelihood function ###
###############################################

LL_observation_term<-function()
{
  return(1)
}

LL_error_term<-function(E,zeta,log=T)
{
  number_of_errors<-0
  number_of_dates<-0
  for(i in 1:length(E))
  {
    number_of_errors<-number_of_errors+sum(E[[i]]==1)
    number_of_dates<-number_of_dates+prod(dim(E[[i]]))
  }
  result<-dbinom(number_of_errors,number_of_dates,zeta,log=log)

  return(result)
}

LL_delays_term<-function()
{
  return(1)
}

