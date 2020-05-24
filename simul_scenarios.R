simul_scenarios <- list()

### Baseline ###
simul_scenarios[["baseline"]] <- list(
  group_size = 100,
  parametric_delays = "gamma",
  ### values for simulation from the NEJM one year [USED]
  # https://www.nejm.org/doi/pdf/10.1056/NEJMc1414992?articleTools=true
  # numbers below are actually from supplement: 
  # https://www.nejm.org/doi/suppl/10.1056/NEJMc1414992/suppl_file/nejmc1414992_appendix.pdf
  mean_delays = data.frame(onset_2_hosp = 5.0,
                           hosp_2_disch = 11.2,
                           hosp_2_death = 4.3,
                           onset_2_report = 5.5),
  sd_delays = data.frame(onset_2_hosp = 4.4,
                         hosp_2_disch = 7.2,
                         hosp_2_death = 4.0,
                         onset_2_report = 5.2),
  prop_missing_data = 0.2,
  prop_error = 0.05,
  error_model = "typo_challenge"
)

### compute CV delays for baseline scenario
simul_scenarios[["baseline"]]$CV_delays <- 
  simul_scenarios[["baseline"]]$sd_delays / 
  simul_scenarios[["baseline"]]$mean_delays

### Low missingness ###

simul_scenarios[["low_missingness"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["low_missingness"]]$prop_missing_data <- 0.05

### Low error ###

simul_scenarios[["low_error"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["low_error"]]$prop_error <- 0.02

### High error ###

simul_scenarios[["high_error"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["high_error"]]$prop_error <- 0.2

### Lognormal delays ###
simul_scenarios[["lognormal_delays"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["lognormal_delays"]]$parametric_delays <- "lognormal"

### Weibull delays ###
simul_scenarios[["weibull_delays"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["weibull_delays"]]$parametric_delays <- "weibull"

### Very small sample size ###
simul_scenarios[["very_small_sample"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["very_small_sample"]]$group_size <- 10

### Small sample size ###
simul_scenarios[["small_sample"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["small_sample"]]$group_size <- 20

### Moderate sample size ###
simul_scenarios[["moderate_sample"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["moderate_sample"]]$group_size <- 50

### Very large sample size ###
simul_scenarios[["very_large_sample"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["very_large_sample"]]$group_size <- 500

### Uniform error ###
simul_scenarios[["uniform_error"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["uniform_error"]]$error_model <- "uniform"

### Long delays ###
simul_scenarios[["long_delays"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["long_delays"]]$mean_delays <- 2 * simul_scenarios[["baseline"]]$mean_delays
# keep CV the same so change sd accordingly
simul_scenarios[["long_delays"]]$sd_delays <- 
  simul_scenarios[["long_delays"]]$mean_delays * 
  simul_scenarios[["long_delays"]]$CV_delays

### Short delays ###
simul_scenarios[["short_delays"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["short_delays"]]$mean_delays <- simul_scenarios[["baseline"]]$mean_delays / 2
# keep CV the same so change sd accordingly
simul_scenarios[["short_delays"]]$sd_delays <- 
  simul_scenarios[["short_delays"]]$mean_delays * 
  simul_scenarios[["short_delays"]]$CV_delays

### High variability in delays ###
simul_scenarios[["highly_variable_delays"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["highly_variable_delays"]]$CV_delays <- 2 * simul_scenarios[["baseline"]]$CV_delays
# change sd accordingly
simul_scenarios[["highly_variable_delays"]]$sd_delays <- 
  simul_scenarios[["highly_variable_delays"]]$mean_delays * 
  simul_scenarios[["highly_variable_delays"]]$CV_delays

### Low variability in delays ###
simul_scenarios[["hardly_variable_delays"]] <- simul_scenarios[["baseline"]]
simul_scenarios[["hardly_variable_delays"]]$CV_delays <- simul_scenarios[["baseline"]]$CV_delays / 2
# change sd accordingly
simul_scenarios[["hardly_variable_delays"]]$sd_delays <- 
  simul_scenarios[["hardly_variable_delays"]]$mean_delays * 
  simul_scenarios[["hardly_variable_delays"]]$CV_delays

### make theta for each scenario

default_index_dates <- function()
{
  index_dates <- list(matrix(c(1, 2), nrow=2), 
                      cbind(c(1, 2), c(1, 3)), 
                      cbind(c(1, 2), c(2, 3), c(1, 4)), 
                      cbind(c(1, 2), c(2, 3), c(1, 4)) )
  names(index_dates) <- c("NoHosp-Alive", "NoHosp-Dead", "Hosp-Alive", "Hosp-Dead")
  index_dates
}

default_index_dates_names <- function()
{
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
  
  index_dates_names
}

set_up_simul <- function(scenario_name, simul_scenarios, n_groups = 4, 
                         n_dates = c(2, 3, 4, 4),
                         range_dates = date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("31/12/2014", "%d/%m/%Y"))),
                         index_dates = default_index_dates(),
                         index_dates_names = default_index_dates_names())
{
  
  theta <- list()
  theta$prop_missing_data <- simul_scenarios[[scenario_name]]$prop_missing_data
  theta$zeta <-  simul_scenarios[[scenario_name]]$prop_error
  theta$mu <- list(simul_scenarios[[scenario_name]]$mean_delays$onset_2_report, 
                   c(simul_scenarios[[scenario_name]]$mean_delays$onset_2_hosp + simul_scenarios[[scenario_name]]$mean_delays$hosp_2_death, simul_scenarios[[scenario_name]]$mean_delays$onset_2_report), 
                   c(simul_scenarios[[scenario_name]]$mean_delays$onset_2_hosp, simul_scenarios[[scenario_name]]$mean_delays$hosp_2_disch, simul_scenarios[[scenario_name]]$mean_delays$onset_2_report), 
                   c(simul_scenarios[[scenario_name]]$mean_delays$onset_2_hosp, simul_scenarios[[scenario_name]]$mean_delays$hosp_2_death, simul_scenarios[[scenario_name]]$mean_delays$onset_2_report))
  theta$CV <- list(simul_scenarios[[scenario_name]]$CV_delays$onset_2_report, 
                   c(simul_scenarios[[scenario_name]]$CV_delays$onset_2_hosp + simul_scenarios[[scenario_name]]$CV_delays$hosp_2_death, simul_scenarios[[scenario_name]]$CV_delays$onset_2_report), 
                   c(simul_scenarios[[scenario_name]]$CV_delays$onset_2_hosp, simul_scenarios[[scenario_name]]$CV_delays$hosp_2_disch, simul_scenarios[[scenario_name]]$CV_delays$onset_2_report), 
                   c(simul_scenarios[[scenario_name]]$CV_delays$onset_2_hosp, simul_scenarios[[scenario_name]]$CV_delays$hosp_2_death, simul_scenarios[[scenario_name]]$CV_delays$onset_2_report))
  n_per_group <- rep(simul_scenarios[[scenario_name]]$group_size, n_groups)
  
  index_dates_order <- compute_index_dates_order(index_dates)
  
  res <- list(theta = theta, n_per_group = n_per_group, range_dates = range_dates, index_dates_names = index_dates_names, delay_dist = simul_scenarios[[scenario_name]]$parametric_delays)
  return(res)
}

### baseline scenario
scenario_name <- "baseline" 
scenario_param <- set_up_simul(scenario_name, simul_scenarios)
# make one simulation 
D <- simul_true_data(scenario_param$theta, scenario_param$n_per_group, scenario_param$range_dates, scenario_param$index_dates_names, delay_dist = scenario_param$delay_dist)
D_with_error <- simul_true_data(scenario_param$theta, scenario_param$n_per_group, scenario_param$range_dates, 
                                scenario_param$index_dates_names, simul_error = TRUE, 
                                remove_allNA_indiv = TRUE, 
                                remove_indiv_at_most_one_date_recorded = TRUE)
tmp <- simul_obs_dat(D$true_dat, scenario_param$theta, scenario_param$range_dates)
obs_dat <- tmp$obs_dat
aug_dat <- list(D=tmp$D, E=tmp$E)
