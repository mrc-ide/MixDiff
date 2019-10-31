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
  name_place_to_load_simulated_data_from <- "83" # "1" # 
  where_to_load_from <- paste0("./SimulatedData/baseline_random_params/",name_place_to_load_simulated_data_from)
  obs_dat <- readRDS(normalizePath(paste0(where_to_load_from,"/SimulatedObsData.rds")))
}

###############################################
### define index_dates ###
###############################################

index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )

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
                       chain_properties=list(n_iter = 5000, burnin = 1, record_every=1),
                        #chain_properties=list(n_iter = 10, burnin = 1, record_every=1),
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
                     index_dates, 
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
plot_parameter_chains(MCMCres, theta_true)
dev.off()

### plot augmented data chains ###
pdf(paste0(where_to_load_from,"/AugDataConvergencePlots_",ext,".pdf"), width=14, height=14)
plot_aug_dat_chains(MCMCres, aug_dat_true)
dev.off()

### get consensus ###

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

posterior_support_for_real_status_of_entry <- function(g, j)
{
  tmp <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][,j] )
  # remove missing entries to focus on true errors, not just missing data
  tmp [tmp == -1] <- NA
  ### definition using the mode posterior -- sesems a bit restrictive 
  #tmp2 <- sapply(seq_len(nrow(MCMCres$aug_dat_chain[[1]]$E[[g]])), function(i) 
  # as.numeric(names(which.max(table(tmp[i,]))) ) == aug_dat_true$E[[g]][i,j] )
  ### definition using threshold 1/4 posterior support
  prob <- sapply(seq_len(nrow(MCMCres$aug_dat_chain[[1]]$E[[g]])), function(i) 
  {
    if(all(is.na(tmp[i,])))
    {
      res <- NA
    } else if (any(tmp[i,] == as.character(aug_dat_true$E[[g]][i,j]))) 
    {
      res <- as.vector(table(tmp[i,])[as.character(aug_dat_true$E[[g]][i,j])]/sum(table(tmp[i,]))) 
    } else res <- 0
    return(res)
  })
  return(prob)
}

find_problematic_Es <- function(g, j, threshold_posterior_support = 1/4)
{
  prob <- which(posterior_support_for_real_status_of_entry(g,j) <= threshold_posterior_support)
  return(prob)
}

posterior_support_list <- lapply(seq_len(length(index_dates)), 
                                 function(g) lapply(seq_len(1+lengths(index_dates)[g]/2), 
                                                    function(j) posterior_support_for_real_status_of_entry(g, j)))

posterior_support <- unlist(posterior_support_list)


par(mfrow = c(2, 1))
hist(posterior_support, col = "grey", breaks = seq(0, 1, 0.01))
hist(posterior_support[posterior_support<0.95], col = "grey", breaks = seq(0, 1, 0.01))

threshold_posterior_support <- 0.25
abline(v = threshold_posterior_support, col = "red", lty = 2)

# using threshold 1/4
prob <- lapply(seq_len(length(index_dates)), 
               function(g) lapply(seq_len(1+lengths(index_dates)[g]/2), 
                                  function(j) find_problematic_Es(g, j, 
                                                                  threshold_posterior_support = threshold_posterior_support)))

# using threshold 1/2
# prob <- lapply(seq_len(length(index_dates)), 
#                function(g) lapply(seq_len(1+lengths(index_dates)[g]/2), 
#                                   function(j) find_problematic_Es(g, j, threshold_posterior_support = 1/2)))
#prob
#prob[[4]]

g <- 3
j <- 3

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
                                     obs_dat, 
                                     posterior = "mode", 
                                     threshold_error_support = 0.95)


par(mfrow = c(2, 2), mar = c(4, 4, 5, .5))
for(g in 1:length(consensus$inferred_E))
{
  image(t(consensus$inferred_E_numeric[[g]]), 
        col = c("black", "forestgreen", "yellow", "orange", "red"), 
        breaks = seq(-0.5, 4.5, 1),
        main = paste("group", g))
}

# save in Excel
write_xlsx_consensus(consensus,
                                 file = "consensus.xlsx",
                                 where = where_to_load_from,
                                 col_width = 10,
                                 overwrite = TRUE)



#####


par(mfrow = c(4, 2))
for(g in 1:4)
{
  image(t(aug_dat_true$E[[g]]), col = c( "red", "forestgreen", "grey"))
  image(t(consensus$consensus_E[[g]]), col = c( "red", "forestgreen", "grey"))
}

g <- 1
plot(aug_dat_true$D[[g]], consensus$D[[g]], col = scales::alpha("grey", 0.5), pch = 19)
for(g in 2:4)
  points(aug_dat_true$D[[g]], consensus$D[[g]], col = scales::alpha("grey", 0.5), pch = 19)
table(aug_dat_true$E[[g]] == consensus$E[[g]])['FALSE'] / sum(table(aug_dat_true$E[[g]] == consensus$E[[g]]))

### Compute sensitivity and specificity of detecting errors in dates
# remove missing dates from the denominator
detec <- MixDiff:::compute_sensitivity_specificity_from_consensus(aug_dat_true, consensus)

### specificity: where we've detected an error which did not exist, what was the reason? 

false_pos <- detec$false_pos
# posterior support for erroneous entry
lapply(1:length(detec$sensitivity), function(g) if(nrow(false_pos[[g]])>0) sapply(1:nrow(false_pos[[g]]), function(i) posterior_support_list[[g]][[false_pos[[g]][i,2]]][false_pos[[g]][i,1]]) else NA )

g <- 1
i <- 18
j <- 1#1
j2 <- 2#2
par(mfrow = c(2, 1))
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, j]), type = "l")
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, j2]), type = "l")

g <- 4
i <- 5
par(mfrow = c(4, 1))
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 1]), type = "l")
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 2]), type = "l")
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 3]), type = "l")
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 4]), type = "l")
aug_dat_true$E[[g]][i,]
aug_dat_true$D[[g]][i,]
obs_dat[[g]][i,]
consensus$E[[g]][i,]

# posterior support for error is not huge
# in group 1 mostly 50% chance of each being an error - which makes sense when only 2 dates observed

### sensitivity: where we've not detected an error, what was the reason? 

false_neg <- detec$false_neg
# posterior support for correct entry
lapply(1:length(detec$sensitivity), function(g) if(nrow(false_neg[[g]])>0) sapply(1:nrow(false_neg[[g]]), function(i) posterior_support_list[[g]][[false_neg[[g]][i,2]]][false_neg[[g]][i,1]]) else NA )

g <- 4
i <- 96#92#85#
j <- 2#4#3#

g <- 3
i <- 88#40#78#
j <- 1#3#3#

g <- 2
i <- 11
j <- 2

g <- 1 # these are the same as the false_pos but for the other date - makes sense
i <- 54
j <- 2
aug_dat_true$D[[g]][i,]
obs_dat[[g]][i,]
aug_dat_true$D[[g]][i,j]
table(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, j]))
hist(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, j]))
aug_dat_true$E[[g]][i,]
consensus$E[[g]][i,]
par(mfrow = c(2, 1))
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, j]), type = "l")
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][i, j]), type = "l")
posterior_support_list[[g]][[false_neg[[g]][which(false_neg[[g]][,1] == i),2]]][false_neg[[g]][which(false_neg[[g]][,1] == i),1]]

g <- 4
i <- 38
ylim <- c(16000, 16300)
par(mfrow = c(1, 1))
plot(sapply(seq_len(length(MCMCres$aug_dat_chain)), 
            function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 1]), 
     type = "l", col = "blue", ylim = ylim)
lines(sapply(seq_len(length(MCMCres$aug_dat_chain)), 
             function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 2]), 
      col = "turquoise")
lines(sapply(seq_len(length(MCMCres$aug_dat_chain)), 
             function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 3]), 
      col = "red")
lines(sapply(seq_len(length(MCMCres$aug_dat_chain)), 
             function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, 4]), 
      col = "purple")
par(xpd = TRUE)
points(rep(length(MCMCres$aug_dat_chain), 4), aug_dat_true$D[[g]][i, ], 
       col = c("blue", "turquoise", "red", "purple"),
       pch = 21, cex = 1.5)
points(rep(length(MCMCres$aug_dat_chain), 4)*1.02, obs_dat[[g]][i, ], 
       col = c("blue", "turquoise", "red", "purple"),
       pch = 19, cex = 1.5)

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

consensus_E_same_as_true_E <- lapply(1:length(consensus$E), function(g) 
{
  apply(consensus$E[[g]] == aug_dat_true$E[[g]], 1, all)
})

id_problematic_consensus <- lapply(consensus_E_same_as_true_E, function(x) which(!x))

### could add something about support - i.e. if consensus is only weakly 
# supported maybe differentiate from cases where consensus is highly supported

### could also identify a priori which individuals are going to be problematic

threshold_consensus_support <- 0.95

g <- 1

one_erroneous_entry <- sapply(1:nrow(aug_dat_true$E[[g]]), function(i) all(sort(aug_dat_true$E[[g]][i,]) == c(0, 1)))
which(one_erroneous_entry)

one_erroneous_entry_one_missing <- sapply(1:nrow(aug_dat_true$E[[g]]), function(i) all(sort(aug_dat_true$E[[g]][i,]) == c(-1, 1)))
which(one_erroneous_entry_one_missing)

two_erroneous_entries <- sapply(1:nrow(aug_dat_true$E[[g]]), function(i) all(sort(aug_dat_true$E[[g]][i,]) == c(1, 1)))
which(two_erroneous_entries)

post_support_this_group <- cbind(posterior_support_list[[1]][[1]], posterior_support_list[[1]][[2]])
low_support_this_group <- post_support_this_group < threshold_consensus_support
tmp <- which(low_support_this_group, arr.ind = TRUE)
tmp[tmp[,2]==1,1]

# so only individual 98 has lower than 0.95 support for one of their dates yet not an a prior problem
aug_dat_true$D[[g]][98,]

g <- 2

post_support_this_group <- cbind(posterior_support_list[[g]][[1]], posterior_support_list[[g]][[2]], posterior_support_list[[g]][[3]])
low_support_this_group <- post_support_this_group < threshold_consensus_support
tmp <- which(low_support_this_group, arr.ind = TRUE)
unique(tmp[,1])

# so all those with problematic consensus have a low posterior support for at least one date

# To measure specificity: are we ever inferring an error for an individual with all observed correct errors
at_least_one_mistake_in_obs <- lapply(1:length(aug_dat_true$E), function(g) 
  sapply(1:nrow(aug_dat_true$E[[g]]), function(i) 
    (1 %in% aug_dat_true$E[[g]][i,])))
at_least_one_mistake_detected_in_consensus <- lapply(1:length(consensus$E), function(g) 
  sapply(1:nrow(consensus$E[[g]]), function(i) 
    (1 %in% consensus$E[[g]][i,])))

true_neg_indiv <- lapply(1:length(aug_dat_true$E), function(g)
  which(at_least_one_mistake_in_obs[[g]] == 0 & at_least_one_mistake_detected_in_consensus[[g]] == 0))

true_pos_indiv <- lapply(1:length(aug_dat_true$E), function(g)
  which(at_least_one_mistake_in_obs[[g]] == 1 & at_least_one_mistake_detected_in_consensus[[g]] == 1))

false_neg_indiv <- lapply(1:length(aug_dat_true$E), function(g)
  which(at_least_one_mistake_in_obs[[g]] == 1 & at_least_one_mistake_detected_in_consensus[[g]] == 0))

false_pos_indiv <- lapply(1:length(aug_dat_true$E), function(g)
  which(at_least_one_mistake_in_obs[[g]] == 0 & at_least_one_mistake_detected_in_consensus[[g]] == 1))

ability_to_detect_individuals_with_errors <- list(
  sensitivity = lengths(true_pos_indiv) / (lengths(true_pos_indiv) + lengths(false_neg_indiv)),
  specifiticy = lengths(true_neg_indiv) / (lengths(true_neg_indiv) + lengths(false_pos_indiv))
)

g <- 4
i <- 5
consensus$E[[g]][i,]
consensus$D[[g]][i,]
aug_dat_true$E[[g]][i,]
aug_dat_true$D[[g]][i,]
obs_dat[[g]][i,]

# issue: sometimes the consensus date will be the true D as here, but the consensus E will be E = 1
# this can happen because for E = 1 there sill be multiple D values, which each will have low support

tmp_chain_E <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][i,1])
plot(tmp_chain_E, typ = "l")
tmp_chain_D <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i,1])
plot(tmp_chain_D, typ = "l")

hist(tmp_chain_E, col = "grey")
hist(tmp_chain_D, col = "grey")

#### examining parameter values ####
are_true_param_in_95perc_post(MCMCres, theta_true)


