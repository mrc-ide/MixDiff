###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

rm(list=ls())

###############################################
### if needed, compile documentation, then check, build and install the package ###
###############################################

roxygen2::roxygenise()
devtools::check()
devtools::build()
devtools::install()

###############################################
### read in data ###
###############################################

USE_SIMULATED_DATA <- TRUE

if(!USE_SIMULATED_DATA)
{
  where_to_load_from <- "./EbolaData"
  raw_dat<-readRDS(paste0(where_to_load_from,"/Dat.rds"))
  colDates <- grep("Date", names(raw_dat))
  tmp <- split(raw_dat[colDates],raw_dat$Path)
  # splitting dataset according to Path and removing NA date columns in each of these
  # - should only remain dates that are relevant for each group
  obs_dat <- lapply(tmp, function(x) sapply(which(colSums(is.na(x))!=nrow(x)), function(j) date_to_int(x[,j]) )) ### converting obs_dat to be integers - easier to handle than dates
} else
{
  name_place_to_load_simulated_data_from <- "1"
  where_to_load_from <- paste0("./SimulatedData/",name_place_to_load_simulated_data_from)
  obs_dat <- readRDS(normalizePath(paste0(where_to_load_from,"/SimulatedObsData.rds")))
}

n_dates <- sapply(obs_dat, ncol )
n_groups <- length(n_dates)

###############################################
### define augmented data to be used for initialisation of the chain ###
###############################################

index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )
aug_dat <- initialise_aug_data(obs_dat, index_dates_order)

###############################################
### define parameters to be used for initialisation of the chain ###
###############################################

index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)

###############################################
### Run the MCMC ###
###############################################

n_iter <- 50000 # currently (21st Nov 2016, updating 1/10th of Di per group at each iteration, 100 iterations take ~360 seconds)

move_D_by_groups_of_size <- 1

### prior parameters 

prior_shape1_prob_error=3
prior_shape2_prob_error=12
prior_mean_mean_delay=100
prior_mean_CV_delay=100

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

logpost_chain <- rep(NA, n_iter)
logpost_chain[1] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay, range_dates)

n_accepted_D_moves <- 0
n_proposed_D_moves <- 0

n_accepted_mu_moves <- 0
n_proposed_mu_moves <- 0

n_accepted_CV_moves <- 0
n_proposed_CV_moves <- 0

#n_accepted_zeta_moves <- 0 # not used as Gibbs sampler 
#n_proposed_zeta_moves <- 0 # not used as Gibbs sampler

### turn on and off various moves, useful for debugging
D_moves_on <- TRUE
mu_moves_on <- TRUE
CV_moves_on <- TRUE
zeta_moves_on <- TRUE

### std of moves

fraction_Di_to_update <- 1/10

sdlog_mu <- 0.15 # for now moving all mus with the same sd, 
# might need to revisit this as some delays might be longer than others an require different sdlog to optimise mixing of the chain

sdlog_CV <- 0.25 # for now moving all CVs with the same sd, 
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
                            prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay, range_dates) 
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
                             prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay) 
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
                                prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay)
          n_proposed_mu_moves <- n_proposed_mu_moves + 1
          n_accepted_mu_moves <- n_accepted_mu_moves + tmp$accept
          if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
        }
      }
    }
    
    # move CV
    if(CV_moves_on)
    {
      for(g in 1:n_groups)
      {
        for(j in 2:ncol(curr_aug_dat$D[[g]]))
        {
          tmp <- move_lognormal(what="CV", g, j-1, sdlog_CV, 
                                curr_aug_dat,
                                curr_theta, 
                                obs_dat, 
                                prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay)
          n_proposed_CV_moves <- n_proposed_CV_moves + 1
          n_accepted_CV_moves <- n_accepted_CV_moves + tmp$accept
          if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
        }
      }
    }
    
    # recording the value of parameters after all moves
    theta_chain <- add_new_value_chain_theta(theta_chain, curr_theta)
    aug_dat_chain <- add_new_value_chain_aug_dat(aug_dat_chain, curr_aug_dat)
    
    # recording the likelihood after all moves
    logpost_chain[k+1] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, prior_shape1_prob_error, prior_shape2_prob_error, prior_mean_mean_delay, prior_mean_CV_delay, range_dates)
  }
})

if(!USE_SIMULATED_DATA)
{
  # add time to name so that can keep track of several results if needed # could change this if want more specific naming of various runs
  ext <- date()
  ext <- gsub(" ","_",ext)
  ext <- gsub(":","",ext)
  save.image(paste0("EbolaData/ResultsEstimation_EbolaData_",ext,".Rdata"))
}else
{
  # add time to name so that can keep track of several results if needed # could change this if want more specific naming of various runs
  ext <- date()
  ext <- gsub(" ","_",ext)
  ext <- gsub(":","",ext)
  save.image(paste0(where_to_load_from,"/ResultsEstimation_SimulatedData_",ext,".Rdata"))
}

###############################################
### acceptance probabilities ###
###############################################

n_accepted_D_moves / n_proposed_D_moves
n_accepted_mu_moves / n_proposed_mu_moves
n_accepted_CV_moves / n_proposed_CV_moves
#n_accepted_zeta_moves / n_proposed_zeta_moves # not computed as now using Gibbs sampler for Zeta

###############################################
### remove burnin ###
###############################################

burnin <- 1:25000
logpost_chain <- logpost_chain[-burnin]
theta_chain$zeta <- theta_chain$zeta[-burnin]
for(g in 1:n_groups)
{
  if(n_dates[g]>=3)
  {
    theta_chain$mu[[g]] <- theta_chain$mu[[g]][-burnin,]
    theta_chain$CV[[g]] <- theta_chain$CV[[g]][-burnin,]
  }else
  {
    theta_chain$mu[[g]] <- theta_chain$mu[[g]][-burnin]
    theta_chain$CV[[g]] <- theta_chain$CV[[g]][-burnin]
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

if(USE_SIMULATED_DATA) theta_simul <- readRDS(paste0(where_to_load_from,"/ThetaUsedForSimulation.rds"))

### parameters ###

pdf(paste0(where_to_load_from,"/ParamConvergencePlots_",ext,".pdf"), width=14, height=7)
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

# looking at CV delay
group_idx <- 1 ##########################
j <- 1
CV <- theta_chain$CV[[group_idx]]
plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(non hospitalised-alive group)", ylim=c(0, 2))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$CV[[group_idx]][j], col=j)
par(xpd=FALSE)
legend("topright", "Onset-Report", lty=1)
group_idx <- 2 ##########################
j <- 1
CV <- theta_chain$CV[[group_idx]][,j]
plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(non hospitalised-dead group)", ylim=c(0, 2))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$CV[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  CV <- theta_chain$CV[[group_idx]][,j]
  lines(CV, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$CV[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 3 ##########################
j <- 1
CV <- theta_chain$CV[[group_idx]][,j]
plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(hospitalised-alive group)", ylim=c(0, 2))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$CV[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  CV <- theta_chain$CV[[group_idx]][,j]
  lines(CV, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$CV[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
group_idx <- 4 ##########################
j <- 1
CV <- theta_chain$CV[[group_idx]][,j]
plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(hospitalised-dead group)", ylim=c(0, 2))
par(xpd=TRUE)
if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$CV[[group_idx]][j], col=j)
par(xpd=FALSE)
for(j in 2:(n_dates[group_idx]-1))
{
  CV <- theta_chain$CV[[group_idx]][,j]
  lines(CV, col=j)
  par(xpd=TRUE)
  if(USE_SIMULATED_DATA) points(n_iter-max(burnin)+n_iter/25,theta_simul$CV[[group_idx]][j], col=j)
  par(xpd=FALSE)
}
legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"), lty=1, col=1:n_dates[group_idx])
dev.off()

### augmented data ###

pdf(paste0(where_to_load_from,"/AugDataConvergencePlots_",ext,".pdf"), width=14, height=14)
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

###############################################
### Correlations ###
###############################################

cor_mu_CV <- list()

par(mfrow=c(2, 5),mar=c(5, 6, 1, 1))
group_idx <- 1
plot(theta_chain$mu[[group_idx]], theta_chain$CV[[group_idx]], type="l")
cor_mu_CV[[group_idx]] <- cor.test(theta_chain$mu[[group_idx]], theta_chain$CV[[group_idx]])

for(group_idx in 2:n_groups)
{
  cor_mu_CV[[group_idx]] <- list()
  for(j in 1:(n_dates[[group_idx]]-1))
  {
    plot(theta_chain$mu[[group_idx]][,j], theta_chain$CV[[group_idx]][,j], type="l", col=j)
    cor_mu_CV[[group_idx]][[j]] <- cor.test(theta_chain$mu[[group_idx]][,j], theta_chain$CV[[group_idx]][,j])
  }
}

cor_mu_CV
# much less correlation now after we have reparameterized to be mean and CV rather than mean and SD

###############################################
### TO DO ###
###############################################

# Anne: 
# check the MCMC, 
# try to speed up if possible
# considering only calculating the likelihood for some iterations (e.g. after burnin and thinning), posthoc? 
# should we update zeta after each D_i move, or after all D_i in a group move? 
# keep track of acceptance rate for D and for mu/CV per group and per deay rather than altogether, to check if some moves are more successful than others. 
# also consider using Gibbs samplers to move mu and CV --> for this need to reformulate as shape/scale: but doesn't seem obvious to sample from the posterior distribution? 
# why do we tend to underestimate the mean delays? related to discretization of gamma distr? 
# write some code to start from last point in the chain
# create a hyperprior list which contains all the prior parameters - easier than keeping track of each of them separately
# in initMCMC.R: index_dates_order A list containing indications on ordering of dates, see details. #### CONSIDER CALCULATING THIS AUTOMATICALLY FROM index_dates

# Marc: 
# finish writing
# think about the 1/(T-T0)

# Future ideas: 
# simulation study
# other datasets - Marc to talk to John? 
# outputs: proportion erroneous data - ... - nice graphs




