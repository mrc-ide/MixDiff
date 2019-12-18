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


