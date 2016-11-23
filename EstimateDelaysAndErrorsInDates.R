###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

rm(list=ls())
library(EpiEstim) # to use DiscrSI which does the discretised Gamma

###############################################
### source functions from other scripts ###
###############################################

source("InitMCMC.R")
source("MCMCMoves.R")
source("Utilities.R")
source("LikelihoodPrior.R")

###############################################
### source objects from Constants.R ###
###############################################

source("Constants.R")

###############################################
### read in data ###
###############################################

USE_SIMULATED_DATA <- TRUE

if(!USE_SIMULATED_DATA)
{
  raw_dat<-readRDS("Dat.rds")
  colDates <- grep("Date", names(raw_dat))
  tmp <- split(raw_dat[colDates],raw_dat$Path)
  # splitting dataset according to Path and removing NA date columns in each of these
  # - should only remain dates that are relevant for each group
  obs_dat <- lapply(tmp, function(x) sapply(which(colSums(is.na(x))!=nrow(x)), function(j) date_to_int(x[,j]) )) ### converting obs_dat to be integers - easier to handle than dates
} else
{
  obs_dat <- readRDS("SimulatedObsData.rds")  
}

n_dates <- sapply(obs_dat, ncol )
n_groups <- length(n_dates)

###############################################
### define parameters to be used for initialisation of the chain ###
###############################################

theta <- initialise_theta(obs_dat, index_dates)

###############################################
### define augmented data to be used for initialisation of the chain ###
###############################################

aug_dat <- initialise_aug_data(obs_dat, index_dates_order)

###############################################
### Run the MCMC ###
###############################################

n_iter <- 1000 # currently (21st Nov 2016, updating 1/10th of Di per group at each iteration, 100 iterations take ~360 seconds)

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

sdlog_mu <- 0.15 # for now moving all mus with the same sd, 
# might need to revisit this as some delays might be longer than others an require different sdlog to optimise mixing of the chain

sdlog_sigma <- 0.25 # for now moving all sigmas with the same sd, 
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

# save.image("ResultsEstimation_EbolaData.Rdata")
# save.image("ResultsEstimation_SimulatedData.Rdata")

###############################################
### acceptance probabilities ###
###############################################

n_accepted_D_moves / n_proposed_D_moves
n_accepted_mu_moves / n_proposed_mu_moves
n_accepted_sigma_moves / n_proposed_sigma_moves

#n_accepted_zeta_moves / n_proposed_zeta_moves # not computed as now using Gibbs sampler for Zeta

###############################################
### remove burnin ###
###############################################

burnin <- 1:500
logpost_chain <- logpost_chain[-burnin]
theta_chain$zeta <- theta_chain$zeta[-burnin]
for(g in 1:n_groups)
{
  if(n_dates[g]>=3)
  {
  theta_chain$mu[[g]] <- theta_chain$mu[[g]][-burnin,]
  theta_chain$sigma[[g]] <- theta_chain$sigma[[g]][-burnin,]
  }else
  {
    theta_chain$mu[[g]] <- theta_chain$mu[[g]][-burnin]
    theta_chain$sigma[[g]] <- theta_chain$sigma[[g]][-burnin]
  }
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

if(USE_SIMULATED_DATA) theta_simul <- readRDS("ThetaUsedForSimulation.rds")

### parameters ###

pdf("ParamConvergencePlots.pdf", width=14, height=7)
par(mfrow=c(2, 5),mar=c(5, 6, 1, 1))

# looking at the logposterior chain 
plot(logpost_chain, type="l", xlab="Iterations", ylab="Log posterior")

# looking at mean delay 
group_idx <- 1 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(non hospitalised-alive group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$mu[[group_idx]][j])
par(xpd=FALSE)

legend("topright", "Onset-Report", lty=1)
group_idx <- 2 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]][,j]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(non hospitalised-dead group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$mu[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  mu <- theta_chain$mu[[group_idx]][,j]
  lines(mu, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$mu[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 3 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]][,j]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(hospitalised-alive group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$mu[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  mu <- theta_chain$mu[[group_idx]][,j]
  lines(mu, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$mu[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 4 ##########################
j <- 1
mu <- theta_chain$mu[[group_idx]][,j]
plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(hospitalised-dead group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$mu[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  mu <- theta_chain$mu[[group_idx]][,j]
  lines(mu, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$mu[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])

# looking at zeta
zeta <- theta_chain$zeta
plot(zeta, type="l", xlab="Iterations", ylab="zeta")
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$zeta)
par(xpd=FALSE)

# looking at std delay
group_idx <- 1 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(non hospitalised-alive group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$sigma[[group_idx]][j], col=j)
par(xpd=FALSE)
legend("topright", "Onset-Report", lty=1)
group_idx <- 2 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]][,j]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(non hospitalised-dead group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$sigma[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  sigma <- theta_chain$sigma[[group_idx]][,j]
  lines(sigma, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$sigma[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 3 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]][,j]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(hospitalised-alive group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$sigma[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  sigma <- theta_chain$sigma[[group_idx]][,j]
  lines(sigma, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$sigma[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 4 ##########################
j <- 1
sigma <- theta_chain$sigma[[group_idx]][,j]
plot(sigma, type="l", xlab="Iterations", ylab="std delays\n(hospitalised-dead group)", ylim=c(0, 20))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$sigma[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  sigma <- theta_chain$sigma[[group_idx]][,j]
  lines(sigma, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$sigma[[group_idx]][j], col=j)
  par(xpd=FALSE)
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
# why do we tend to underestimate the mean delays? related to discretization of gamma distr? 
# write some code to start from last point in the chain
# check correlations between outputs

# Marc: 
# finish writing
# think about the 1/(T-T0)

# Future ideas: 
# simulation study
# other datasets - Marc to talk to John? 
# outputs: proportion erroneous data - ... - nice graphs




