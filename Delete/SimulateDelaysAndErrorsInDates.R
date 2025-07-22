rm(list=ls())

###############################################
### if needed, compile documentation, then check, build and install the package ###
###############################################

devtools::document()
devtools::check()
devtools::build()
devtools::install()

####################################
### creating a simulated dataset ###
####################################

n_groups <- 4
n_dates <- c(2, 3, 4, 4)

theta <- list()
theta$prop_missing_data <- 0.05 ### this is missing from the estimation model (directly available from the data so not explicitely modelled)
theta$zeta <- 0.05 ### probability that, when not missing, the date is recorded with error
theta$mu <- list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13))
theta$CV <- list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))
n_per_group <- rep(100, n_groups)
range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
index_dates_order <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 3), c(1, 4)) )

D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
D_with_error <- simul_true_data(theta, n_per_group, range_dates, index_dates, simul_error = TRUE, remove_allNA_indiv=TRUE)
tmp <- simul_obs_dat(D$true_dat, theta, range_dates)
E <- tmp$E
obs_dat <- tmp$obs_dat

aug_dat <- list(D=D$true_dat, E=tmp$E)

####################################
### saving this somewhere ###
####################################

name_place_to_save <- "2" # in the future find a way to automatically generate a random folder name

where_to_save <- paste0("./SimulatedData/",name_place_to_save)
if(!dir.exists("./SimulatedData/")) dir.create("./SimulatedData/")
if(!dir.exists(where_to_save)) dir.create(where_to_save)

saveRDS(obs_dat, file = normalizePath(paste0(where_to_save,"/SimulatedObsData.rds")))
saveRDS(aug_dat, file = normalizePath(paste0(where_to_save,"/SimulatedAugData.rds")))
saveRDS(theta, file = normalizePath(paste0(where_to_save,"/ThetaUsedForSimulation.rds")))



