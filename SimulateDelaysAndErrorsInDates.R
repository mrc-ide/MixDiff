rm(list=ls())

###############################################
### if needed, compile documentation, then check, build and install the package ###
###############################################

roxygen2::roxygenise()
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

####################################
### THIS IS WHERE I AM IN CHECKING STUFF WORK - DOCUMENTING FUNCTIONS APPROPRIATELY ###
####################################

tmp <- simul_obs_dat(D, theta, range_dates)
E <- tmp$E
obs_dat <- tmp$obs_dat
# remove those with only missing dates - assuming that you always have at least one date present
for(g in 1:n_groups)
{
  exclude <- which(rowSums(is.na(obs_dat[[g]]))==ncol(obs_dat[[g]]))
  if(length(exclude)>0)
  {
    obs_dat[[g]] <- obs_dat[[g]][-exclude,]
    D[[g]] <- D[[g]][-exclude,]
    E[[g]] <- E[[g]][-exclude,]
  }
}
aug_dat <- list(D=D, E=E)
rm(D)
rm(E)

####################################
### saving this somewhere ###
####################################

name_place_to_save <- "1" # in the future find a way to automatically generate a random folder name

where_to_save <- paste0("./SimulatedData/",name_place_to_save)
if(!dir.exists("./SimulatedData/")) dir.create("./SimulatedData/")
if(!dir.exists(where_to_save)) dir.create(where_to_save)

saveRDS(obs_dat, file = normalizePath(paste0(where_to_save,"/SimulatedObsData.rds")))
saveRDS(aug_dat, file = normalizePath(paste0(where_to_save,"/SimulatedAugData.rds")))
saveRDS(theta, file = normalizePath(paste0(where_to_save,"/ThetaUsedForSimulation.rds")))



