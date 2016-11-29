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

###############################################
### define index_dates andindex_dates_order  ###
###############################################

index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
### NOTE THE FOLLOWING SHOULD BE AUTOMATICALLY CALCULATED FROM index_dates ###
index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )

###############################################
### MCMC settings ###
###############################################

MCMC_settings <- list( moves_switch=list(D_on = TRUE, mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
                       moves_options=list(fraction_Di_to_update = 1/10, move_D_by_groups_of_size = 1, sdlog_mu = 0.15, sdlog_CV = 0.25), 
                       chain_properties=list(n_iter = 100, burnin = 50, record_every=10))
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

MCMCres <- RunMCMC(obs_dat, 
                   MCMC_settings,
                   hyperpriors,
                   index_dates,
                   index_dates_order) ### CHANGE THIS SO index_dates_order is computed automatically from index_dates

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
###############################################
### THIS IS WHERE I AM AT IN TERMS OF RECODING INTO THE PACKAGE
###############################################
###############################################



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
# in initMCMC.R: index_dates_order A list containing indications on ordering of dates, see details. #### CONSIDER CALCULATING THIS AUTOMATICALLY FROM index_dates
# currently initialisation of augmented data can still start in stupid place, where order of dates is ok but delays are too large, so make sure this is not the case
# everywhere replace 1:n by seq_len(n) 
# where we use ncol(curr_aug_dat$D[[g]]), check this as I think it may need to be defined from index_dates rather than from D

# Marc: 
# finish writing
# think about the 1/(T-T0)

# Future ideas: 
# simulation study
# other datasets - Marc to talk to John? 
# outputs: proportion erroneous data - ... - nice graphs




