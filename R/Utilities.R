
###############################################
### functions to add a new value at end of chains ###
###############################################

add_new_value_chain_theta <- function(theta_chain, new_theta)
{
  n_groups <- length(theta_chain$mu)
  theta_chain$zeta <- c(theta_chain$zeta, new_theta$zeta)
  for(g in 1:n_groups)
  {
    theta_chain$mu[[g]] <- rbind(theta_chain$mu[[g]], new_theta$mu[[g]])
    theta_chain$CV[[g]] <- rbind(theta_chain$CV[[g]], new_theta$CV[[g]])
  }
  return(theta_chain)
}

add_new_value_chain_aug_dat <- function(curr_aug_dat, new_aug_dat)
{
  n_groups <- length(curr_aug_dat$D)
  for(g in 1:n_groups)
  {
    for(j in 1:length(curr_aug_dat$D[[g]]))
    {
      curr_aug_dat$D[[g]][[j]] <- rbind(curr_aug_dat$D[[g]][[j]], new_aug_dat$D[[g]][,j])
      curr_aug_dat$E[[g]][[j]] <- rbind(curr_aug_dat$E[[g]][[j]], new_aug_dat$E[[g]][,j])
    }
  }
  return(curr_aug_dat)
}

###############################################
### functions to handle dates ###
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
### functions to find parameters of distributions from mean/var or mean/sd ###
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

find_params_gamma <- function(mean, sigma) # function to determine parameters of the gamma distribution corresponding to a given mean and std
{
  shape <- (mean/sigma)^2
  scale <- sigma^2/(mean)
  return(c(shape, scale))
}
# test if works: 
# param_gamma <- find_params_gamma(0.1, 0.05)
# sample_gamma <- rgamma(1000, shape=param_gamma[1], scale=param_gamma[2])
# mean(sample_gamma)
# sd(sample_gamma)

find_params_gamma_from_mean_CV <- function(mean, CV) # function to determine parameters of the gamma distribution corresponding to a given mean and CV
{
  sigma <- mean*CV
  res <- find_params_gamma(mean, sigma)
  return(res)
}
# test if works: 
# param_gamma <- find_params_gamma_from_mean_CV(0.1, 0.5)
# sample_gamma <- rgamma(1000, shape=param_gamma[1], scale=param_gamma[2])
# mean(sample_gamma)
# sd(sample_gamma)/mean(sample_gamma)

###############################################
### compute_delta functions to compute relevant delays based on index, which tells you which dates should be used for delayl calculation ###
###############################################

compute_delta_group_delay_and_indiv<-function(D, group_idx, delay_idx, indiv_idx, index = index_dates)
{
  Delta <- D[[group_idx]][indiv_idx,index[[group_idx]][,delay_idx][2]] - D[[group_idx]][indiv_idx,index[[group_idx]][,delay_idx][1]]
  return(Delta)
}

compute_delta <- function(D, index = index_dates)
{
  Delta <- lapply(1:length(D), function(g){
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
### find range of dates
###############################################

find_range <- function(obs_dat)
{
  min_date <- min(obs_dat[[1]][,1], na.rm=TRUE)
  max_date <- max(obs_dat[[1]][,1], na.rm=TRUE)
  for(g in 1:length(obs_dat))
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