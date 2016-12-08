###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

rm(list=ls())

###############################################
### if needed, compile documentation, then check, build and install the package ###
###############################################

devtools::document()
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

###############################################
### define index_dates andindex_dates_order  ###
###############################################

index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
### NOTE THE FOLLOWING SHOULD BE AUTOMATICALLY CALCULATED FROM index_dates ###
index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )

###############################################
### MCMC settings ###
###############################################

MCMC_settings <- list( moves_switch=list(D_on = TRUE, E_on = TRUE,  swapE_on = TRUE,  mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
                       moves_options=list(fraction_Di_to_update = 1/10, move_D_by_groups_of_size = 1, fraction_Ei_to_update = 1/10, sdlog_mu = 0.15, sdlog_CV = 0.25), 
                       init_options=list(mindelay=0, maxdelay=100),
                       #chain_properties=list(n_iter = 200, burnin = 1, record_every=1))
                       chain_properties=list(n_iter = 1000, burnin = 250, record_every=2))
#chain_properties=list(n_iter = 5000, burnin = 500, record_every=10))
#chain_properties=list(n_iter = 50000, burnin = 5000, record_every=50))
#chain_properties=list(n_iter = 250000, burnin = 50000, record_every=100))
# for now moving all mus and CVs with the same sd, 
# might need to revisit this as some delays might be longer than others an require different sdlog to optimise mixing of the chain

###############################################
### prior parameters  ###
###############################################

hyperpriors <- list(
  shape1_prob_error=3, 
  shape2_prob_error=12, 
  mean_mean_delay=100, 
  mean_CV_delay=100)

###############################################
### Run the MCMC  ###
###############################################

set.seed(1)
#Rprof()
system.time({
  #profvis::profvis({
  MCMCres <- RunMCMC(obs_dat, 
                     MCMC_settings,
                     hyperpriors,
                     index_dates,
                     index_dates_order) ### CHANGE THIS SO index_dates_order is computed automatically from index_dates
})
#Rprof(NULL)
#summaryRprof()
# 2 Dec --> n_iter = 5000, burnin = 500, record_every=10 takes 596secs


###############################################
### save results ###
###############################################

if(!USE_SIMULATED_DATA)
{
  # add time to name so that can keep track of several results if needed # could change this if want more specific naming of various runs
  ext <- date()
  ext <- gsub(" ","_",ext)
  ext <- gsub(":","",ext)
  saveRDS(MCMCres, paste0("EbolaData/ResultsEstimation_EbolaData_",ext,".rds"))
}else
{
  # add time to name so that can keep track of several results if needed # could change this if want more specific naming of various runs
  ext <- date()
  ext <- gsub(" ","_",ext)
  ext <- gsub(":","",ext)
  saveRDS(MCMCres, paste0(where_to_load_from,"/ResultsEstimation_SimulatedData_",ext,".rds"))
}

###############################################
### plotting the MCMC output ###
###############################################

# If working on simulated data, load the parameters used for simulation for comparison with MCMC estimates
if(USE_SIMULATED_DATA) 
{
  theta_true <- readRDS(paste0(where_to_load_from,"/ThetaUsedForSimulation.rds"))
  aug_dat_true <- readRDS(paste0(where_to_load_from,"/SimulatedAugData.rds"))
}else
{
  theta_true <- NULL
  aug_dat_true <- NULL
}

### plot parameter chains ###
pdf(paste0(where_to_load_from,"/ParamConvergencePlots_",ext,".pdf"), width=14, height=7)
plot_parameter_chains(MCMCres, theta_true)
dev.off()

### plot augmented data chains ###
pdf(paste0(where_to_load_from,"/AugDataConvergencePlots_",ext,".pdf"), width=14, height=14)
plot_aug_dat_chains(MCMCres, aug_dat_true)
dev.off()

###############################################
### Correlations and autocorrelation ###
###############################################

cor_mu_CV <- compute_correlations_mu_CV(MCMCres)
autocorr <- compute_autocorr(MCMCres)

###############################################
### Get and plot posterior estimates ###
###############################################

pdf(paste0(where_to_load_from,"/PosteriorDistrPlots_",ext,".pdf"), width=14, height=7)
MCMCres_summary <- get_param_posterior_estimates(MCMCres, theta_true=theta_true, cex.axis=0.8)
dev.off()

###############################################
### Investigating issue with estimation of mean delay from onset to hosp in group 3 ###
###############################################

### checking the actual 'true' simulated data

par(mfrow=c(1,3))
hist(aug_dat_true$D[[3]][,2] - aug_dat_true$D[[3]][,1])
summary(aug_dat_true$D[[3]][,2] - aug_dat_true$D[[3]][,1])
hist(aug_dat_true$D[[3]][,3] - aug_dat_true$D[[3]][,2])
summary(aug_dat_true$D[[3]][,3] - aug_dat_true$D[[3]][,2])
hist(aug_dat_true$D[[3]][,4] - aug_dat_true$D[[3]][,1])
summary(aug_dat_true$D[[3]][,4] - aug_dat_true$D[[3]][,1])

### comparing the posterior distribution of true aug_dat and parameters to posterior reached by MCMC chain

# posterior of true aug data and true parameters
lposterior_total(aug_dat_true, theta_true, obs_dat, hyperpriors, index_dates, range_dates=NULL)

# posterior of current aug data and true parameters
lposterior_total(MCMCres$aug_dat_chain[[length(MCMCres$aug_dat_chain)]], theta_true, obs_dat, hyperpriors, index_dates, range_dates=NULL)

# highest posterior reached by MCMC chain
max(MCMCres$logpost_chain)

###############################################
### Examining how well we reestimate the E (error/missingness in data) ###
###############################################

par(mfrow=c(4, 2))
dy <- 10
for(g in 1:4)
{
  tmp <- t(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) table(MCMCres$aug_dat_chain[[e]]$E[[g]])))
  ylim <- mean(tmp[,2])+c(-1,1)*dy
  plot(tmp[,2], type="l", main="Number of correctly recorded dates in that group", ylab="", xlab="Iterations", ylim=ylim)
  abline(h=table(aug_dat_true$E[[g]])[2], col="red")
  ylim <- mean(tmp[,3])+c(-1,1)*dy
  plot(tmp[,3], type="l", main="Number of erroroneous dates in that group", ylab="", xlab="Iterations", ylim=ylim)
  abline(h=table(aug_dat_true$E[[g]])[3], col="red")
}

### define inferred E as wrong if support for true E < 1/4 

find_problematic_Es <- function(g, j)
{
  tmp <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][,j] )
  ### definition using the mode posterior -- sesems a bit restrictive 
  #tmp2 <- sapply(seq_len(nrow(MCMCres$aug_dat_chain[[1]]$E[[g]])), function(i) as.numeric(names(which.max(table(tmp[i,]))) ) == aug_dat_true$E[[g]][i,j] )
  ### definition using threshold 1/4 posterior support
  threshold_posterior_support <- 1/4
  tmp2 <- sapply(seq_len(nrow(MCMCres$aug_dat_chain[[1]]$E[[g]])), function(i) as.vector(table(tmp[i,])[as.character(aug_dat_true$E[[g]][i,j])]/sum(table(tmp[i,]))>threshold_posterior_support ))
  prob <- which(!tmp2)
  return(prob)
}

prob <- lapply(seq_len(length(index_dates)), function(g) lapply(seq_len(1+lengths(index_dates)[g]/2), function(j) find_problematic_Es(g, j)))

prob
prob[[4]]

g <- 4
#j <- 2
j <- 4


prob_i <- prob[[g]][[j]][1]
sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][prob_i, j])
table(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][prob_i, j]))
table(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][prob_i, j]))
hist(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][prob_i, j]))
obs_dat[[g]][prob_i, j]
obs_dat[[g]][prob_i, ]
aug_dat_true$E[[g]][prob_i,]
aug_dat_true$D[[g]][prob_i,]
MCMCres$aug_dat_chain[[length(MCMCres$aug_dat_chain)]]$E[[g]][prob_i,]
MCMCres$aug_dat_chain[[length(MCMCres$aug_dat_chain)]]$D[[g]][prob_i,]


### All of the remaining issues seem to be because the observation is erroneous but the error is too small to be detected, we can't do anything about this.

###############################################
### TO DO ###
###############################################

# Anne: 
# check the MCMC, 
# try to speed up if possible
# keep track of acceptance rate for D and for mu/CV per group and per delay rather than altogether, to check if some moves are more successful than others. 
# write some code to start from last point in the chain
# in initMCMC.R: index_dates_order A list containing indications on ordering of dates, see details. #### CONSIDER CALCULATING THIS AUTOMATICALLY FROM index_dates
# where we use ncol(curr_aug_dat$D[[g]]), check this as I think it may need to be defined from index_dates rather than from D
# question for Rich: should all functions used in tests be "public"?
# do we indeed want to update zeta after each D_i move? maybe not useful? 

# Marc: 
# finish writing
# think about the 1/(T-T0)

# Future ideas: 
# simulation study
# other datasets - Marc to talk to John? 
# outputs: proportion erroneous data - ... - nice graphs




