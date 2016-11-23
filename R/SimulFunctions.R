#######################################
### functions to simulate a dataset ###
#######################################

simul_true_data <- function(theta, n_per_group, range_dates, index_dates)
{
  D <- list() 
  for(g in 1:length(theta$mu))
  {
    D[[g]] <- matrix(NA, n_per_group[g], length(theta$mu[[g]])+1)
    D[[g]][,1] <- sample(range_dates[1]:range_dates[2], n_per_group[g], replace = TRUE)
    for(j in 1:ncol(index_dates[[g]]))
    {
      params <- find_params_gamma_from_mean_CV(theta$mu[[g]][j], theta$CV[[g]][j])
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


