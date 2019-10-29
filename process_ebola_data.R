library(MixDiff)

### read in data from Tini's paper
### NOTE: one of the issues in this dataset is that I am not sure the date of outcome recorded is actually the date of death or the date of discharge...

where_to_load_from <- "./EbolaData/WestAfrica_Garske"

dat <- read.csv(paste0(where_to_load_from, "/data.csv"), 
                header = TRUE, stringsAsFactors = FALSE)

### add identifier

dat$id <- 1:nrow(dat)

### convert dates to dates:

dat_cols <- grep("Date", names(dat))
names(dat)[dat_cols]
for(kk in dat_cols)
{
  dat[,kk] <- as.Date(dat[,kk])
}

### generate a variable for our four groups depending on final status and hospitalisation status

dat$group_name <- paste0(dat$HospitalizedEver, dat$FinalStatus)
dat$group_name[grep("NA", dat$group_name)] <- NA

dat$group <- NA
dat$group[dat$group_name %in% "NoAlive"] <- 1
dat$group[dat$group_name %in% "NoDead"] <- 2
dat$group[dat$group_name %in% "YesAlive"] <- 3
dat$group[dat$group_name %in% "YesDead"] <- 4

### choose a way of subsetting the data to get reasonable numbers for analysis to run in decent amount of time

dat_with_unknown_groups <- dat
dat <- dat[!is.na(dat$group),]

table(dat$QuarterOnsetInferred, dat$Country)
table(dat$group, dat$Country)

dat_Liberia_early2014 <- dat[
  dat$QuarterOnsetInferred %in% c("Jan - Mar 2014", "Apr - Jun 2014", "Jul - Sep 2014") & dat$Country %in% "Liberia",]
table(dat_Liberia_early2014$group)

### format the data as needed to apply MixDiff

dat_chosen <- dat_Liberia_early2014

obs_dat <- list(dat_chosen[dat_chosen$group %in% 1, c("DateOnset", "DateReport")],
                dat_chosen[dat_chosen$group %in% 2, c("DateOnset", "DateOutcomeComp", "DateReport")],
                dat_chosen[dat_chosen$group %in% 3, c("DateOnset", "DateHospitalCurrentAdmit", "DateOutcomeComp", "DateReport")],
                dat_chosen[dat_chosen$group %in% 4, c("DateOnset", "DateHospitalCurrentAdmit", "DateOutcomeComp", "DateReport")])

id_obs_dat <- list(dat_chosen[dat_chosen$group %in% 1, "id"],
                   dat_chosen[dat_chosen$group %in% 2, "id"],
                   dat_chosen[dat_chosen$group %in% 3, "id"],
                   dat_chosen[dat_chosen$group %in% 4, "id"])

obs_dat <- lapply(obs_dat, function(x) {
  sapply(1:ncol(x), function(e) MixDiff::date_to_int(x[,e]))
})

### remove any individual with no dates recorded

all_NAs <- lapply(obs_dat, function(x) {
  which(sapply(1:nrow(x), function(e) all(is.na(x[e,]))))
})
obs_dat <- lapply(1:length(obs_dat), function(g) {
  res <- obs_dat[[g]]
  if(any(all_NAs[[g]])) res <- res[-all_NAs[[g]],]
  return(res)
})
id_obs_dat <- lapply(1:length(id_obs_dat), function(g) {
  res <- id_obs_dat[[g]]
  if(any(all_NAs[[g]])) res <- res[-all_NAs[[g]]]
  return(res)
})

sapply(obs_dat, dim) 

### take only the first 100 of each group

obs_dat_large <- obs_dat
id_obs_dat_large <- id_obs_dat
obs_dat <- lapply(obs_dat, function(x) x[1:100,])
id_obs_dat <- lapply(id_obs_dat, function(x) x[1:100])

sapply(obs_dat, dim) 

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
                       #chain_properties=list(n_iter = 10, burnin = 1, record_every=1),
                       #chain_properties=list(n_iter = 500, burnin = 1, record_every=1),
                       chain_properties=list(n_iter = 5000, burnin = 500, record_every=10),
                       tol = 1e-6)

###############################################
### prior parameters  ###
###############################################

hyperparameters <- list(
  shape1_prob_error=3, 
  shape2_prob_error=12, 
  mean_mean_delay=100, 
  mean_CV_delay=100) # AT THE MOMENT NOT USED 

MCMCres <- RunMCMC(obs_dat, 
                   MCMC_settings,
                   hyperparameters,
                   index_dates, 
                   seed = 1)

ext <- date()
ext <- gsub(" ","_",ext)
ext <- gsub(":","",ext)
ext2 <- paste0(ext, "_",
               MCMC_settings$chain_properties$n_iter, "iter_",
               MCMC_settings$chain_properties$burnin, "burnt_",
               MCMC_settings$chain_properties$record_every, "thin_")
saveRDS(MCMCres, paste0("./EbolaData/WestAfrica_Garske/ResultsEstimation_SimulatedData_",ext2,".rds"))

### plot parameter chains ###
pdf(paste0(where_to_load_from,"/ParamConvergencePlots_",ext,".pdf"), width=14, height=7)
plot_parameter_chains(MCMCres)
dev.off()

### plot augmented data chains ###
pdf(paste0(where_to_load_from,"/AugDataConvergencePlots_",ext,".pdf"), width=14, height=14)
plot_aug_dat_chains(MCMCres)
dev.off()

### Get and plot posterior estimates ###
pdf(paste0(where_to_load_from,"/PosteriorDistrPlots_",ext,".pdf"), width=14, height=7)
MCMCres_summary <- get_param_posterior_estimates(MCMCres, cex.axis=0.8)
dev.off()

### Get and plot consensus ###
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
                     file = "/consensus.xlsx",
                     where = where_to_load_from,
                     col_width = 10,
                     overwrite = TRUE)


### look at the errors we found ###
### for this would be good to keep the IDs of cases so we could check more stuff in initial dataset
lapply(1:length(consensus$inferred_E), function(g) table(consensus$inferred_E[[g]]))

g <- 3
idx <- which(consensus$inferred_E[[g]] == "error_high_support", arr.ind = TRUE)
ii <- idx[1,1]
# what was observed for that case
obs_dat[[g]][ii, ]
int_to_date(obs_dat[[g]][ii, ])
# what was inferred for that case
consensus$inferred_D[[g]][ii,]
int_to_date(consensus$inferred_D[[g]][ii,])
# original dataset:
dat[dat$id %in% id_obs_dat[[g]][ii],]
# check convergence: 
tmp_D <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][ii, ])
tmp_E <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][ii, ])
par(mfrow = c(2, nrow(tmp_E)))
for(k in 1:nrow(tmp_E))
  plot(tmp_E[k,], type = "l", ylim = c(-1, 1), col = k)
for(k in 1:nrow(tmp_E))
  plot(tmp_D[k,], type = "l", ylim = range(tmp_D), col = k)

g <- 4
idx <- which(consensus$inferred_E[[g]] == "error_high_support", arr.ind = TRUE)
ii <- idx[1,1]
# what was observed for that case
obs_dat[[g]][ii, ]
int_to_date(obs_dat[[g]][ii, ])
# what was inferred for that case
consensus$inferred_D[[g]][ii,]
int_to_date(consensus$inferred_D[[g]][ii,])
# original dataset:
dat[dat$id %in% id_obs_dat[[g]][ii],]
# check convergence: 
tmp_D <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][ii, ])
tmp_E <- sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][ii, ])
par(mfrow = c(2, nrow(tmp_E)))
for(k in 1:nrow(tmp_E))
  plot(tmp_E[k,], type = "l", ylim = c(-1, 1), col = k)
for(k in 1:nrow(tmp_E))
  plot(tmp_D[k,], type = "l", ylim = range(tmp_D), col = k)

### compare the onset we infer to what was inferred in Tini's dataset ###

# for those with missing onset date

manually_inferred <- NULL
automatically_inferred <- NULL
ids_inferred <- NULL

for(g in 1:length(obs_dat))
{
  missing_onset <- is.na(obs_dat[[g]][,1])
  if(any(missing_onset))
  {
    ids_missing_onset <- id_obs_dat[[g]][which(missing_onset)] 
    ids_inferred <- c(ids_inferred, ids_missing_onset)
    manually_inferred <- c(manually_inferred, dat$DateOnsetInferred[dat$id %in% ids_missing_onset])
    automatically_inferred <- c(automatically_inferred, int_to_date(consensus$inferred_D[[g]][which(missing_onset)] ))
  }
}

#differences between automatic inference and manual inference: 
diff_inferred <- automatically_inferred - manually_inferred 
names(diff_inferred) <- ids_inferred
# check the very big ones
dat[dat$id %in% names(diff_inferred)[1],]
dat[dat$id %in% names(diff_inferred)[3],]

par(mfrow = c(1, 1))
plot(manually_inferred, automatically_inferred)
abline(0, 1)

### seems like we estimate a much shorter onset to report delay in group 2 
### so this may explain the differences







