rm(list=ls())

###############################################
### if needed, compile documentation, then check, build and install the package ###
###############################################

devtools::document()
devtools::check()
devtools::build()
devtools::install()

####################################
### creating a set of dummy simulated datasets ###
####################################

# n_groups <- 4
# n_dates <- c(2, 3, 4, 4)
# 
# theta <- list()
# theta$prop_missing_data <- 0.05 ### this is missing from the estimation model (directly available from the data so not explicitely modelled)
# theta$zeta <- 0.05 ### probability that, when not missing, the date is recorded with error
# theta$mu <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
# theta$CV <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))
# n_per_group <- rep(100, n_groups)
# range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
# index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
# index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )
# 
# for(name_place_to_save in 1:100)
# {
#   
#   D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#   D_with_error <- simul_true_data(theta, n_per_group, range_dates, index_dates, simul_error = TRUE, remove_allNA_indiv=TRUE)
#   tmp <- simul_obs_dat(D$true_dat, theta, range_dates)
#   E <- tmp$E
#   obs_dat <- tmp$obs_dat
#   
#   aug_dat <- list(D=D$true_dat, E=tmp$E)
#   
#   ####################################
#   ### saving this somewhere ###
#   ####################################
#   
#   where_to_save <- paste0("./SimulatedData/baseline_random_params/",name_place_to_save)
#   if(!dir.exists("./SimulatedData/")) dir.create("./SimulatedData/baseline_random_params/")
#   if(!dir.exists(where_to_save)) dir.create(where_to_save)
#   
#   saveRDS(obs_dat, file = normalizePath(paste0(where_to_save,"/SimulatedObsData.rds")))
#   saveRDS(aug_dat, file = normalizePath(paste0(where_to_save,"/SimulatedAugData.rds")))
#   saveRDS(theta, file = normalizePath(paste0(where_to_save,"/ThetaUsedForSimulation.rds")))
#   
# }



####################################
### creating a set of ebola like simulated datasets ###
####################################

### values for simulation from NEJM 9 months [NOT USED HERE]
# https://www.nejm.org/doi/pdf/10.1056/NEJMoa1411100?articleTools=true
# mean / sd for all countries

### onset to hosp 
# 5.0 / 4.7

### hosp to discharge
# 11.8 / 6.1

### hosp to death
# 4.2 / 6.4

### onset to report
# 6.1 / 8.5


### values for simulation from the NEJM one year [USED]
# https://www.nejm.org/doi/pdf/10.1056/NEJMc1414992?articleTools=true
# numbers below are actually from supplement: 
# https://www.nejm.org/doi/suppl/10.1056/NEJMc1414992/suppl_file/nejmc1414992_appendix.pdf

# mean / CI (for the mean) for all countries, Dec 2013 to 25 Nov 2014
# sd / CI (for the sd) for all countries, Dec 2013 to 25 Nov 2014
# for confirmed and probable cases

### onset to hosp 
# 5.0 / 4.9-5.1
mean_onset_2_hosp <- 5.0
# 4.4 / 4.3-4.6
sd_onset_2_hosp <- 4.4

### hosp to discharge
# 11.2 / 10.8-11.7
mean_hosp_2_disch <- 11.2
# 7.2 / 6.8-7.6
sd_hosp_2_disch <- 7.2

### hosp to death
# 4.3 / 4.1-4.5
mean_hosp_2_death <- 4.3
# 4.0 / 3.8-4.3
sd_hosp_2_death <- 4.0

### onset to report (not reported in main text?)
# 5.5 / 5.4-5.7
mean_onset_2_report <- 5.5
# 5.2 / 5.1-5.3
sd_onset_2_report <- 5.2

n_groups <- 4
n_dates <- c(2, 3, 4, 4)

theta <- list()
theta$prop_missing_data <- 0.2 ### this is the proportion of missing data - in Tini's Ebola in West Africa Dataset this is 22% 
theta$zeta <- 0.05 ### probability that, when not missing, the date is recorded with error --> taken roughly as in typo challenge proportion of dates entered wrongly but which are still dates
theta$mu <- list(mean_onset_2_report, 
                 c(mean_onset_2_hosp + mean_hosp_2_death, mean_onset_2_report), 
                 c(mean_onset_2_hosp, mean_hosp_2_disch, mean_onset_2_report), 
                 c(mean_onset_2_hosp, mean_hosp_2_death, mean_onset_2_report))
theta$CV <- list(sd_onset_2_report/mean_onset_2_report, 
                 c(sqrt(sd_onset_2_hosp^2 + sd_hosp_2_death^2)/(mean_onset_2_hosp + mean_hosp_2_death), sd_onset_2_report/mean_onset_2_report), 
                 c(sd_onset_2_hosp/mean_onset_2_hosp, sd_hosp_2_disch/mean_hosp_2_disch, sd_onset_2_report/mean_onset_2_report), 
                 c(sd_onset_2_hosp/mean_onset_2_hosp, sd_hosp_2_death/mean_hosp_2_death, sd_onset_2_report/mean_onset_2_report))
n_per_group <- rep(100, n_groups)
range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("31/12/2014", "%d/%m/%Y")))

index_dates <- list(matrix(c(1, 2), nrow=2), 
                    cbind(c(1, 2), c(1, 3)), 
                    cbind(c(1, 2), c(2, 3), c(1, 4)), 
                    cbind(c(1, 2), c(2, 3), c(1, 4)) )
names(index_dates) <- c("NoHosp-Alive", "NoHosp-Dead", "Hosp-Alive", "Hosp-Dead")

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

index_dates_order <- compute_index_dates_order(index_dates)
  
for(name_place_to_save in 1:100)
{
  
  D <- simul_true_data(theta, n_per_group, range_dates, index_dates_names)
  # D_with_error <- simul_true_data(theta, n_per_group, range_dates, 
  #                                 index_dates_names, simul_error = TRUE, 
  #                                 remove_allNA_indiv=TRUE)
  D_with_error <- simul_true_data(theta, n_per_group, range_dates, 
                                  index_dates_names, simul_error = TRUE, 
                                  remove_allNA_indiv = TRUE, 
                                  remove_indiv_at_most_one_date_recorded=TRUE)
  tmp <- simul_obs_dat(D$true_dat, theta, range_dates)
  obs_dat <- tmp$obs_dat
  aug_dat <- list(D=tmp$D, E=tmp$E)
  
  ####################################
  ### saving this somewhere ###
  ####################################
  
  where_to_save <- paste0("./SimulatedData/baseline_ebola_like/",name_place_to_save)
  if(!dir.exists("./SimulatedData/baseline_ebola_like/")) dir.create("./SimulatedData/baseline_ebola_like/")
  if(!dir.exists(where_to_save)) dir.create(where_to_save)
  
  saveRDS(obs_dat, file = normalizePath(paste0(where_to_save,"/SimulatedObsData.rds")))
  saveRDS(aug_dat, file = normalizePath(paste0(where_to_save,"/SimulatedAugData.rds")))
  saveRDS(theta, file = normalizePath(paste0(where_to_save,"/ThetaUsedForSimulation.rds")))
  
}



