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