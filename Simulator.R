
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

#######################################
### functions to simulate a dataset ###
#######################################

find_params_gamma <- function(mu, sigma)
{
  a <- ((mu - 1)/sigma)^2
  b <- sigma^2/(mu - 1)
  return(c(a, b))
}

simul_true_data <- function(theta, n_per_group, range_dates, index_dates)
{
  D <- list()
  for(g in 1:length(theta$mu))
  {
    D[[g]] <- matrix(NA, n_per_group[g], length(theta$mu[[g]])+1)
    D[[g]][,1] <- sample(range_dates[1]:range_dates[2], n_per_group[g], replace = TRUE)
    for(j in 1:ncol(index_dates[[g]]))
    {
      params <- find_params_gamma(theta$mu[[g]], theta$sigma[[g]])
      delay <- rgamma(n_per_group[g], shape=params[1], scale=params[2])
      D[[g]][,index_dates[[g]][2,j]]  <- D[[g]][,index_dates[[g]][1,j]] + round(delay)
    }
  }
  return(D)
}

simul_obs_dat <- function(D, theta)
{
  E <- D
  obs_dat <- D
  for(g in 1:length(D))
  {
    for(j in 1:ncol(D[[g]]))
    {
      E[[g]][,j] <- sample(c(-1, 1, 0), nrow(D[[g]]), replace=TRUE, prob=c(theta$prop_missing_data, (1-theta$prop_missing_data)*theta$zeta, (1-theta$prop_missing_data)*(1-theta$zeta)))
      obs_dat[[g]][E[[g]][,j]==-1,j]  <- NA
      obs_dat[[g]][E[[g]][,j]==0,j]  <- D[[g]][E[[g]][,j]==0,j]
      obs_dat[[g]][E[[g]][,j]==1,j]  <- sample(range_dates[1]:range_dates[2], sum(E[[g]][,j]==1), replace = TRUE) # need to update if change error model
    }
  }
  return(list(obs_dat=obs_dat, E=E))
}

####################################
### creating a simulated dataset ###
####################################

n_groups <- 4
n_dates <- c(2, 3, 4, 4)

theta <- list()
theta$prop_missing_data <- 0.1 ### this is currently missing from the estimation model
theta$zeta <- 0.1 ### probability that, when not missing, the date is recorded with error
theta$mu <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
theta$sigma <- list(3, c(3, 3), c(3, 3, 3), c(3, 3, 3))
n_per_group <- rep(100, 4)
range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )


D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
tmp <- simul_obs_dat(D, theta)
E <- tmp$E
obs_dat <- tmp$obs_dat
aug_dat <- list(D=D, E=E)
saveRDS(obs_dat, file = "SimulatedObsData.rds")
saveRDS(aug_dat, file = "SimulatedAugData.rds")

#########################################
### applying the MCMC to this dataset ###
#########################################

# execute code from EstimateDelaysAndErrorsInDates.R

# ISSUE obs_data is integers not dates