devtools::load_all()

sim_data <- readRDS("./simulated_baseline/sim_data.rds")
sim_estim_baseline <- readRDS("./simulated_baseline/sim_estim_baseline.rds")

g <- 4 # fourth group
index_dates <- sim_estim_baseline$baseline[[1]]$index_dates[[g]]

n_simul <- length(sim_data$baseline)

mean_true_delay <- mean_obs_delay <- matrix(NA, n_simul, ncol(index_dates))

for(i in 1:n_simul) {
  #i <- 1 # first simulation

  true_delays <- compute_delta(D = sim_data$baseline[[i]]$true_dat,
                               index_dates = sim_estim_baseline$baseline[[i]]$index_dates)
  mean_true_delay[i, ] <- apply(true_delays[[g]], 2, mean)

  obs_delays <- compute_delta(D = sim_data$baseline[[i]]$obs_dat,
                              index_dates = sim_estim_baseline$baseline[[i]]$index_dates)
  mean_obs_delay[i, ] <- apply(obs_delays[[g]], 2, mean, na.rm = TRUE)
}

colMeans(mean_true_delay)
colMeans(mean_obs_delay)

summary(mean_true_delay)
summary(mean_obs_delay)

###

plot(sim_estim_baseline$baseline[[i]]$logpost_chain, type = "l")

mu_chain <- lapply(sim_estim_baseline$baseline, function(ee) t(sapply(ee$theta_chain, function(e) e$mu[[g]])))

par(mfrow = c(1, ncol(index_dates)))

for(d in 1:ncol(index_dates)) {
  plot(mu_chain[[1]][, d], type = "l", col = scales::alpha("black", .1), ylim = c(0, 20))
  for(i in 2:n_simul) {
    lines(mu_chain[[i]][, d], type = "l", col = scales::alpha("black", .1))
  }
  abline(h = c(11, 12, 13)[d], col = "red", lty = 2)
}

d_est <- t(sapply(mu_chain, function(e) apply(e, 2, mean)))
for(d in 1:ncol(index_dates)) {
  hist(d_est[, d], breaks = seq(0, 20, .1))
  abline(v = c(11, 12, 13)[d], col = "red", lty = 2)
}

coverage_est <- sapply(1:ncol(index_dates), function(d) sum(sapply(mu_chain, function(e) quantile(e[, d], 0.025) < c(11, 12, 13)[d] & c(11, 12, 13)[d] < quantile(e[, d], 0.975))) / n_simul)
coverage_est

## example of problematic one

i <- 4

D1_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$D[[g]][, 1]))
D2_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$D[[g]][, 2]))
D3_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$D[[g]][, 3]))
D4_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$D[[g]][, 4]))
E1_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$E[[g]][, 1]))
E2_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$E[[g]][, 2]))
E3_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$E[[g]][, 3]))
E4_chain <- t(sapply(sim_estim_baseline$baseline[[i]]$aug_dat_chain, function(e) e$E[[g]][, 4]))

D1_true <- sim_data$baseline[[i]]$true_dat[[g]][, 1]
D2_true <- sim_data$baseline[[i]]$true_dat[[g]][, 2]
D3_true <- sim_data$baseline[[i]]$true_dat[[g]][, 3]
D4_true <- sim_data$baseline[[i]]$true_dat[[g]][, 4]

par(mfrow = c(5, 5), mar = c(1.5, .5, 1.5, .5))
for(j in 1:25) {
  plot(D1_chain[, j], type = "l", lwd = 2)
  abline(h = D1_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}
for(j in 26:50) {
  plot(D1_chain[, j], type = "l", lwd = 2)
  abline(h = D1_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}
for(j in 51:75) {
  plot(D1_chain[, j], type = "l", lwd = 2)
  abline(h = D1_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}
for(j in 76:100) {
  plot(D1_chain[, j], type = "l", lwd = 2)
  abline(h = D1_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}

par(mfrow = c(5, 5), mar = c(1.5, .5, 1.5, .5))
for(j in 1:25) {
  plot(D2_chain[, j], type = "l", lwd = 2)
  abline(h = D2_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}
for(j in 26:50) {
  plot(D2_chain[, j], type = "l", lwd = 2)
  abline(h = D2_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}
for(j in 51:75) {
  plot(D2_chain[, j], type = "l", lwd = 2)
  abline(h = D2_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}
for(j in 76:100) {
  plot(D2_chain[, j], type = "l", lwd = 2)
  abline(h = D2_true[j], col = "red", lty = 2, lwd = 2)
  title(main = j)
}

# D1 wrong
j = 9
j = 14
j = 53
j = 66
j = 97

# D2 wrong
j = 9
j = 69
j = 97

# Investigate those

j = 9 # swap_E didn't seem to have worked
j = 53 # same
j = 69
sim_data$baseline[[i]]$true_dat[[g]][j, ]
sim_data$baseline[[i]]$obs_dat[[g]][j, ]

par(mfrow = c(2, 2), mar = c(3, 3, 3, 3))
plot(D1_chain[, j], type = "l", lwd = 2, ylim = range(c(D1_chain[, j], D1_true[j])))
abline(h = D1_true[j], col = "red", lty = 2, lwd = 2)
E1_chain[, j]
plot(D2_chain[, j], type = "l", lwd = 2, ylim = range(c(D2_chain[, j], D2_true[j])))
abline(h = D2_true[j], col = "red", lty = 2, lwd = 2)
E2_chain[, j]
plot(D3_chain[, j], type = "l", lwd = 2, ylim = range(c(D3_chain[, j], D3_true[j])))
abline(h = D3_true[j], col = "red", lty = 2, lwd = 2)
E3_chain[, j]
plot(D4_chain[, j], type = "l", lwd = 2, ylim = range(c(D4_chain[, j], D4_true[j])))
abline(h = D4_true[j], col = "red", lty = 2, lwd = 2)
E4_chain[, j]

## these are examples where the swap doesn't look like it has worked: j = 75, j = 93

## what I would expect for j = 66
sim_data$baseline[[i]]$obs_dat[[g]][j, 2] - 11
sim_data$baseline[[i]]$obs_dat[[g]][j, 4] - 13

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## get things ready to put into swap_E
curr_aug_dat = sim_estim_baseline$baseline[[4]]$aug_dat_chain[[1]]
theta = sim_estim_baseline$baseline[[4]]$theta_chain[[1]]
obs_dat = sim_data$baseline[[4]]$obs_dat
hyperparameters <- list(
  # scalars giving the 1st and 2nd shape parameters for the beta prior for zeta
  shape1_prob_error = 3,
  shape2_prob_error = 12,
  # scalars giving the mean of the exponential prior used for mu and CV
  mean_mean_delay = 100,
  mean_CV_delay = 100
)
range_dates = NULL
i <- 9 # this is one of the examples of swap not working
group_idx = g
index_dates <- sim_estim_baseline$baseline[[1]]$index_dates
set.seed(1)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Problems currently when two wrongs one right
obs_dat[[group_idx]][i,]
curr_aug_dat$D[[group_idx]][i,]
curr_aug_dat$E[[group_idx]][i,]
proposed_aug_dat_intermediate$D[[group_idx]][i,]
proposed_aug_dat_intermediate$E[[group_idx]][i,]
proposed_aug_dat$D[[group_idx]][i,]
proposed_aug_dat$E[[group_idx]][i,]
proposed_aug_dat_fin$D[[group_idx]][i,]
proposed_aug_dat_fin$E[[group_idx]][i,]
proposed_aug_dat_rev1$D[[group_idx]][i,]
proposed_aug_dat_rev1$E[[group_idx]][i,]
proposed_aug_dat_rev2$D[[group_idx]][i,]
proposed_aug_dat_rev2$E[[group_idx]][i,]
proposed_aug_dat_rev3$D[[group_idx]][i,]
proposed_aug_dat_rev3$E[[group_idx]][i,]
##

# MCMC_settings <- list(
#   moves_switch = list(D_on = TRUE, E_on = TRUE, swapE_on = TRUE,
#                       mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
#   moves_options = list(
#     fraction_Di_to_update = 1 / 10,
#     move_D_by_groups_of_size = 1,
#     fraction_Ei_to_update = 1 / 10,
#     sdlog_mu = list(
#       0.05,
#       c(0.15, 0.15),
#       c(0.15, 0.15, 0.15),
#       c(0.25, 0.25, 0.25)
#     ),
#     sdlog_CV = list(
#       0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
#   ),
#   init_options = list(
#     mindelay = 0,
#     maxdelay = 20
#   ),
#   chain_properties = list(
#     n_iter = 500,
#     burnin = 50,
#     record_every = 10
#   )
# )
