###############################################
###############################################
# TO DO:
# apply to Ebola data BEFORE cleaning dates, i.e check that format is fine but not that order etc makes sense before applying.
# also: idea to apply to only the early cases as we are not accounting for changes in distributions over time (e.g. prompter hospitalisation later on)
###############################################
###############################################

###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

rm(list=ls())

###############################################
### if needed, compile documentation, then check, build and install the package ###
###############################################

### to install from local version (ANNE)
##devtools::install_github("mrc-ide/rmnist")
#devtools::document()
#devtools::check()
#devtools::build()
#devtools::install()

### to install from remote (NATSUKO and others)
# devtools::install_github("MJomaba/MixDiff@anne")

### load library
library(MixDiff)

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
  #where_to_load_from <- paste0("./SimulatedData/baseline_random_params/",name_place_to_load_simulated_data_from)
  where_to_load_from <- paste0("./SimulatedData/baseline_ebola_like/",name_place_to_load_simulated_data_from)
  obs_dat <- readRDS(normalizePath(paste0(where_to_load_from,"/SimulatedObsData.rds")))
  
  ### add group names and column names ### in the future best to do this in simul function
  group_names <- c("NoHosp-Alive", "NoHosp-Dead", "Hosp-Alive", "Hosp-Dead")
  names(obs_dat) <- group_names
  
  colnames(obs_dat[[1]]) <- c("onset", "report")
  colnames(obs_dat[[2]]) <- c("onset", "death", "report")
  colnames(obs_dat[[3]]) <- c("onset", "hosp.", "discharge", "report")
  colnames(obs_dat[[4]]) <- c("onset", "hosp.", "death", "report")
  
  ### add ids for cases ### in the future best to do this in simul function
  for(g in 1:length(obs_dat))
  {
    rownames(obs_dat[[g]]) <- paste(group_names[[g]], 1:nrow(obs_dat[[g]]), sep = "_")
  }
  
}

###############################################
### define index_dates ###
###############################################

index_dates <- list(matrix(c(1, 2), nrow=2), 
                    cbind(c(1, 2), c(1, 3)), 
                    cbind(c(1, 2), c(2, 3), c(1, 4)), 
                    cbind(c(1, 2), c(2, 3), c(1, 4)) )

names(index_dates) <- names(obs_dat)

index_dates_names <- index_dates
index_dates_names[["NoHosp-Alive"]][,1] <- c("onset", "report")
index_dates_names[["NoHosp-Dead"]][,1] <- c("onset", "death")
index_dates_names[["NoHosp-Dead"]][,2] <- c("onset", "report")
index_dates_names[["Hosp-Alive"]][,1] <- c("onset", "hosp.")
index_dates_names[["Hosp-Alive"]][,2] <- c("hosp.", "discharge")
index_dates_names[["Hosp-Alive"]][,3] <- c("onset", "report")
index_dates_names[["Hosp-Dead"]][,1] <- c("onset", "hosp.")
index_dates_names[["Hosp-Dead"]][,2] <- c("hosp.", "death")
index_dates_names[["Hosp-Dead"]][,3] <- c("onset", "report")

for(g in 1:length(index_dates_names))
  colnames(index_dates_names[[g]]) <- paste("delay", 1:ncol(index_dates_names[[g]]), sep = "_")

###############################################
### MCMC settings ###
###############################################

MCMC_settings <- list( moves_switch=list(D_on = TRUE, E_on = TRUE,  swapE_on = TRUE,  mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
                       moves_options=list(fraction_Di_to_update = 1/10, move_D_by_groups_of_size = 1, fraction_Ei_to_update = 1/10, 
                                          sdlog_mu = list(0.05, c(0.15, 0.15), c(0.15, 0.15, 0.15), c(0.25, 0.25, 0.25)), 
                                          sdlog_CV = list(0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))), 
                       init_options=list(mindelay=0, maxdelay=100),
                       #                       chain_properties=list(n_iter = 200, burnin = 1, record_every=1))
                       #chain_properties=list(n_iter = 500, burnin = 1, record_every=1),
                       #chain_properties=list(n_iter = 500, burnin = 250, record_every=2),
                       #chain_properties=list(n_iter = 5000, burnin = 500, record_every=10),
                       #chain_properties=list(n_iter = 5000, burnin = 1, record_every=1),
                       chain_properties=list(n_iter = 10, burnin = 1, record_every=1),
                       tol = 1e-6)
#chain_properties=list(n_iter = 50000, burnin = 5000, record_every=50))
#chain_properties=list(n_iter = 250000, burnin = 50000, record_every=100))
# for now moving all mus and CVs with the same sd, 
# might need to revisit this as some delays might be longer than others an require different sdlog to optimise mixing of the chain

###############################################
### prior parameters  ###
###############################################

hyperparameters <- list(
  shape1_prob_error=3, 
  shape2_prob_error=12, 
  mean_mean_delay=100, 
  mean_CV_delay=100) # AT THE MOMENT NOT USED 

###############################################
### Run the MCMC  ###
###############################################

#Rprof()
system.time({
  #  profvis::profvis({
  MCMCres <- RunMCMC(obs_dat, 
                     MCMC_settings,
                     hyperparameters,
                     index_dates = index_dates_names, 
                     seed = 2)
})
#Rprof(NULL)
#summaryRprof()
# 2 Dec --> n_iter = 5000, burnin = 500, record_every=10 takes 596secs
# 21 Oct --> n_iter = 5000, burnin = 500, record_every=10 takes 35mins

###############################################
### save results ###
###############################################

if(!USE_SIMULATED_DATA)
{
  # add time to name so that can keep track of several results if needed # could change this if want more specific naming of various runs
  ext <- date()
  ext <- gsub(" ","_",ext)
  ext <- gsub(":","",ext)
  ext2 <- paste0(ext, "_",
                 MCMC_settings$chain_properties$n_iter, "iter_",
                 MCMC_settings$chain_properties$burnin, "burnt_",
                 MCMC_settings$chain_properties$record_every, "thin_")
  saveRDS(MCMCres, paste0("EbolaData/ResultsEstimation_EbolaData_",ext2,".rds"))
}else
{
  # add time to name so that can keep track of several results if needed # could change this if want more specific naming of various runs
  ext <- date()
  ext <- gsub(" ","_",ext)
  ext <- gsub(":","",ext)
  ext2 <- paste0(ext, "_",
                 MCMC_settings$chain_properties$n_iter, "iter_",
                 MCMC_settings$chain_properties$burnin, "burnt_",
                 MCMC_settings$chain_properties$record_every, "thin_")
  saveRDS(MCMCres, paste0(where_to_load_from,"/ResultsEstimation_SimulatedData_",ext2,".rds"))
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
plot_parameter_chains(MCMCres, theta_true, index_dates_names)
dev.off()

### plot augmented data chains ###
pdf(paste0(where_to_load_from,"/AugDataConvergencePlots_",ext,".pdf"), width=14, height=14)
plot_aug_dat_chains(MCMCres, aug_dat_true)
dev.off()

### get consensus ###

###############################################
### Correlations and autocorrelation ###
###############################################

cor_mu_CV <- compute_correlations_mu_CV(MCMCres, index_dates = index_dates_names)
autocorr <- compute_autocorr(MCMCres, index_dates = index_dates_names)

###############################################
### Get and plot posterior estimates ###
###############################################

pdf(paste0(where_to_load_from,"/PosteriorDistrPlots_",ext,".pdf"), width=14, height=7)
MCMCres_summary <- get_param_posterior_estimates(MCMCres, 
                                                 index_dates = index_dates_names, 
                                                 theta_true=theta_true, 
                                                 cex.axis=0.8)
dev.off()

###############################################
### Investigating issue with estimation of mean delay from onset to hosp in group 3 ###
###############################################

### checking the actual 'true' simulated data

par(mfrow=c(1,3), mar = c(5, 5, 5, 1))
hist(aug_dat_true$D[[3]][,2] - aug_dat_true$D[[3]][,1], 
     col = "grey", main = "True delays 1 (Date 1 to Date 2), group 3")
summary(aug_dat_true$D[[3]][,2] - aug_dat_true$D[[3]][,1])
hist(aug_dat_true$D[[3]][,3] - aug_dat_true$D[[3]][,2],
     col = "grey", main = "True delays 2 (Date 2 to Date 3), group 3")
summary(aug_dat_true$D[[3]][,3] - aug_dat_true$D[[3]][,2])
hist(aug_dat_true$D[[3]][,4] - aug_dat_true$D[[3]][,1],
     col = "grey", main = "True delays 3 (Date 1 to Date 4), group 3")
summary(aug_dat_true$D[[3]][,4] - aug_dat_true$D[[3]][,1])

### comparing the posterior distribution of true aug_dat and parameters to posterior reached by MCMC chain

# posterior of true aug data and true parameters
lposterior_total(aug_dat_true, theta_true, obs_dat, hyperparameters, index_dates, range_dates=NULL)

# posterior of current aug data and true parameters
lposterior_total(MCMCres$aug_dat_chain[[length(MCMCres$aug_dat_chain)]], theta_true, obs_dat, hyperparameters, index_dates, range_dates=NULL)
#lposterior_total(MCMCres$aug_dat_chain[[length(MCMCres$aug_dat_chain)]], MCMCres$theta_chain[[length(MCMCres$theta_chain)]], obs_dat, hyperparameters, index_dates, range_dates=NULL)

# highest posterior reached by MCMC chain
max(MCMCres$logpost_chain)

###############################################
### TO DO ###
###############################################

# Anne: 
# check the MCMC, [DO WE RECOVER PRIOR IF WE REMOVE LIKELIHOOD?]
# try to speed up if possible
# keep track of acceptance rate for D and for mu/CV per group and per delay rather than altogether, to check if some moves are more successful than others. 
# write some code to start from last point in the chain
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

consensus <- get_consensus(MCMCres, 
                           posterior = "mode")

inferred <- get_inferred_from_consensus(consensus, 
                                        threshold_error_support = 0.95)

thresholds <- seq(0.5, 1, 0.05)
inferred_all_thresholds <- lapply(thresholds, function(t) 
  get_inferred_from_consensus(consensus, 
                              threshold_error_support = t)) 
names(inferred_all_thresholds) <- thresholds

par(mfrow = c(2, 2), mar = c(4, 4, 5, .5))
for(g in 1:length(inferred$inferred_E))
{
  image(t(inferred$inferred_E_numeric[[g]]), 
        col = c("black", "forestgreen", "yellow", "orange", "red"), 
        breaks = seq(-0.5, 4.5, 1),
        main = paste("group", g))
}

# save in Excel
write_xlsx_inferred(inferred_all_thresholds[["0.95"]],
                    file = "consensus.xlsx",
                    where = where_to_load_from,
                    col_width = 10,
                    overwrite = TRUE)

# visualise
par(mfrow = c(4, 2), mar = c(3, 3, 1, .5))
for(g in 1:4)
{
  image(t(aug_dat_true$E[[g]]), col = c("grey", "white", "red"))
  image(t(inferred_all_thresholds[["0.95"]]$inferred_E_numeric[[g]]), col = c("grey", "white", "yellow", "orange", "red"))
}

# visualise the support for the erroneous / non erroneous dates according to consensus
# can also help decide on thresholds in the future
plot_support_for_error(aug_dat_true, consensus)

### Compute sensitivity and specificity of detecting errors in dates
# remove missing dates from the denominator
detec_dates <- compute_performance_per_date_from_inferred(aug_dat_true, inferred_all_thresholds[["0.95"]])

detec_dates_all_thresholds <- lapply(inferred_all_thresholds, function(e) 
  compute_performance_per_date_from_inferred(aug_dat_true, e))
names(detec_dates_all_thresholds) <- thresholds

sensitivity_dates_all_thresholds <- sapply(detec_dates_all_thresholds, function(e) 
  e$sensitivity)
specificity_dates_all_thresholds <- sapply(detec_dates_all_thresholds, function(e) 
  e$specificity)

ROC_dates <- ROC_per_date(MCMCres, aug_dat_true, thresholds)

plot_ROC <- function(ROC_dates, 
                     mfrow = c(2, 2), 
                     xlim = c(0, 1), ylim = c(0, 1), 
                     ...) # other graphical parameters for plot()
{
  if(nrow(ROC_dates$sensitivity) > prod(mfrow))
  {
    warning("All groups cannot be plotted, adjust mfrow parameters to show all groups")
  }
  par(mfrow = mfrow)
  for(g in 1:nrow(ROC_dates$sensitivity))
  {
    plot(1 - ROC_dates$specificity[g,], ROC_dates$sensitivity[g,], 
         type = "o",
         xlim = xlim, ylim = ylim,
         ...)
  }
}

plot_ROC(ROC_dates, xlim = c(0, 0.1), ylim = c(0, 1))

###### IDEAS: 
# ROC curve to show how threshold for posterior support affects 
# sensitivity / specificity???

###### TO DO: 
# check whether we need some correction factor on the proposal with the new 
# swap function which also changes the missing data
# could improve move_Di by drawing directly from the mixture of delays, as we now do in swap E_i
# add correction factor for move_Di? I think we need one 

###### Remaining issues:
### 1. When there are only two dates recorded and one is erroneous, 
# it's impossisble to know which one of the two is wrong
# so based on the consensus date 
# we will sometimes detect an error where there is none (lack of specificity)
# (this is the ONLY instance I can see of lack of specificity)
# and we will sometimes NOT detect an error which exists (lack of sensitivity)
# HOWEVER in these cases the posterior support for the error status is low, 
# usually around 50% as we oscillate between error for one or the other date

### 2. When the only recorded date is erroneous (because other ones are missing), 
# it's impossible to know that the recorded date is incorrect
# so we loose sensitivity here

### 3. When the erroneous date is too close to the true date, 
# we also can't detect it
# so again loss of sensitivity

### I think in the example above 
# ("Wed_Oct__9_110650_2019_5000iter_500burnt_10thin_")
# this covers ALL reasons for imperfect performance

###############################################################
# Looking at individuals rather than dates
###############################################################

detec_indiv <- compute_performance_per_individual_from_inferred(aug_dat_true, inferred_all_thresholds[["0.95"]])

detec_indiv_all_thresholds <- lapply(inferred_all_thresholds, function(e) 
  compute_performance_per_individual_from_inferred(aug_dat_true, e))
names(detec_indiv_all_thresholds) <- thresholds

sensitivity_indiv_all_thresholds <- sapply(detec_indiv_all_thresholds, function(e) 
  e$sensitivity)
specificity_indiv_all_thresholds <- sapply(detec_indiv_all_thresholds, function(e) 
  e$specificity)

ROC_indiv <- ROC_per_individual(MCMCres, aug_dat_true, thresholds)

plot_ROC(ROC_dates, xlim = c(0, 0.1), ylim = c(0, 1))
plot_ROC(ROC_indiv, xlim = c(0, 0.1), ylim = c(0, 1))

#### examining parameter values ####
are_true_param_in_95perc_post(MCMCres, theta_true)


