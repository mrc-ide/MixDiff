#' Runs the MCMC estimation procedure.
#'
#' @param obs_dat A list of observed data, containing dates `D` and error
#'  indicators `E` for each group.
#' @param MCMC_settings A list of settings to be used for running the MCMC, see
#'  details.
#' @param hyperparameters A list of hyperparameters for zeta, mu and CV. See
#'  details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#'
#' @details \code{MCMC_settings} should be a list containing:
#' \itemize{
#'  \item{\code{moves_switch} : A list of booleans (D_on ,E_on, mu_on, CV_on,
#'   zeta_on) stating whether each parameter/augmented data should be moved in
#'    the procedure or not.}
#'  \item{\code{moves_options}: A list of the following elements:
#'  \itemize{
#'  \item{\code{fraction_Di_to_update}: The fraction of augmented dates to be
#'   updated at each iteration of the MCMC.}
#'  \item{\code{move_D_by_groups_of_size}: The number of augmented dates to be
#'   updated simultaneously in each group.}
#'  \item{\code{fraction_Ei_to_update}: The fraction of indicators of whether
#'   observed dates are erroneous to be updated at each iteration of the MCMC.}
#'  \item{\code{sdlog_mu}: The standard deviations to be used for proposing
#'   moves of the mean delays. This should be a list of length
#'    \code{n_groups = length(obs_dat)}. Each element in the list should be a
#'     vector with length given by the numbers of delays to be considered in
#'      this group.}
#'  \item{\code{sdlog_CV}: The standard deviations to be used for proposing
#'   moves of the CV of delays. This should be a list of length
#'    \code{n_groups = length(obs_dat)}. Each element in the list should be a
#'     vector with length given by the numbers of delays to be considered in
#'      this group.}
#'  }
#'  }
#'  \item{\code{init_options}: A list of the following elements:
#'  \itemize{
#'  \item{\code{mindelay}: The minimum delay, below which dates are considered
#'   incompatile with one another at the initialisation stage of the MCMC.}
#'  \item{\code{maxdelay}: The maximum delay, above which dates are considered
#'   incompatile with one another at the initialisation stage of the MCMC.}
#'  \item{\code{record_every}: A number indicating, after the burnin, every
#'   how many iterations outputs should be recorded.}
#'   }
#'  }
#'  \item{\code{chain_properties}: A list of the following elements:
#'  \itemize{
#'  \item{\code{n_iter}: The total number of iteration of MCMC to run.}
#'  \item{\code{burnin}: The number of initial iterations to consider as the
#'   burnin period - no output is recorded for these initial MCMC iterations.}
#'  \item{\code{record_every}: A number indicating, after the burnin, every
#'   how many iterations outputs should be recorded.}
#'   }
#'  }
#' }
#' \code{hyperparameters} should be a list containing:
#' \itemize{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape parameter
#'   for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential
#'   prior used for parameter \code{theta$CV}}
#' }
#' \code{index_dates} should be a list of length
#'  \code{n_groups = length(obs_dat)}. Each element of \code{index_dates} should
#'   be a matrix with 2 rows and a number of columns corresponding to the delays
#'    of interest for that group. For each column (i.e. each delay), the first
#'     row gives the index of the origin date, and the second row gives the
#'      index of the destination date.
#' The number of columns of index_dates[[k]] should match the length of
#'  theta$mu[[k]] and theta$CV[[k]]
#'
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and
#'  c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are
#'   respectively the mean and coefficient of variation of two delays: the
#'    first delay being between date 1 and date 2, and the second being between
#'     date 1 and date 3.
#'
#' @return A list of the following elements:
#'  \itemize{
#'  \item{\code{theta_chain}: a list of parameters (mu, cv, zeta), at each
#'   recorded step of the MCMC chain}
#'  \item{\code{aug_dat_chain}: a list of augmented dates and error indicators,
#'   at each recorded step of the MCMC chain}
#'  \item{\code{logpost_chain}: a vector of values of the log posterior at
#'   each recorded step of the MCMC chain}
#'  \item{\code{accept_prob}: a list of the probabilities of acceptance for
#'   each move type across all MCMC iterations}
#' }
#'
#' @export
#'
#' @examples
#' # Simulate data to use
#' n_groups <- 4
#' index_dates <- list(
#'   matrix(c(1, 2), nrow = 2),
#'   cbind(c(1, 2), c(1, 3)),
#'   cbind(c(1, 2), c(2, 3),c(1, 4)),
#'   cbind(c(1, 2), c(2, 3), c(1, 4))
#' )
#' theta <- list(
#'   mu = list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13)),
#'   CV = list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5)),
#'   prop_missing_data = 0.2,
#'   zeta = 0.05
#'  )
#'
#' n_per_group <- rep(10, n_groups)
#' range_dates <- c(0, 30)
#'
#' simul_dat <- simul_true_data(theta, n_per_group, range_dates, index_dates,
#'                              simul_error = TRUE)
#' obs_dat <- simul_dat$obs_dat
#'
#' # Set up hyperparameters
#' hyperparameters <- list(
#'     shape1_prob_error = 3,
#'     shape2_prob_error = 12,
#'     mean_mean_delay = 10,
#'     mean_CV_delay = 10)
#'
#' # Set up MCMC
#' MCMC_settings <- list(
#'   moves_switch = list(D_on = TRUE, E_on = TRUE, swapE_on = TRUE,
#'                     mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
#'   moves_options = list(
#'     fraction_Di_to_update = 1 / 10,
#'     move_D_by_groups_of_size = 1,
#'     fraction_Ei_to_update = 1 / 10,
#'     sdlog_mu = list(
#'       0.05,
#'       c(0.15, 0.15),
#'       c(0.15, 0.15, 0.15),
#'       c(0.25, 0.25, 0.25)
#'     ),
#'     sdlog_CV = list(
#'       0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
#'   ),
#'   init_options = list(
#'     mindelay = 0,
#'     maxdelay = 20
#'   ),
#'   chain_properties = list(
#'     n_iter = 500,
#'     burnin = 50,
#'     record_every = 10
#'   )
#' )
#'
#' # Run MCMC
#' MCMC_result <- RunMCMC(obs_dat,
#'                        MCMC_settings,
#'                        hyperparameters,
#'                        index_dates)
#'
#' # Results
#' MCMC_result$theta_chain
#' MCMC_result$aug_dat_chain
#' MCMC_result$logpost_chain
#' plot(MCMC_result$logpost_chain, type = "l",
#'      ylab = "Log Posterior", xlab = "Iteration")
#' MCMC_result$accept_prob
#'
RunMCMC <- function(obs_dat,
                    MCMC_settings,
                    hyperparameters,
                    index_dates) {

  # Initialise dimensions and check inputs
  n_dates <- sapply(obs_dat, ncol)
  n_groups <- length(n_dates)

  ## Note: need to add checks for obs_dat input e.g. check for all NA obs_dat rows.
  ## Make sure error messages are informative.

  check_MCMC_settings(MCMC_settings, index_dates)

  #----------------------------------------------------------------------------
  # Initialise augmented data and parameters

  aug_dat <- initialise_aug_data(obs_dat,
                                 index_dates, #compute_index_dates_order(index_dates),
                                 MCMC_settings)
  theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)
  range_dates <- find_range(obs_dat)

  #----------------------------------------------------------------------------
  # Initialise storage of MCMC chains

  curr_theta <- theta
  theta_chain <- list()
  theta_chain[[1]] <- curr_theta

  curr_aug_dat <- aug_dat
  aug_dat_chain <- list()
  aug_dat_chain[[1]] <- curr_aug_dat

  logpost_chain <- rep(NA, (MCMC_settings$chain_properties$n_iter -
                            MCMC_settings$chain_properties$burnin) /
                            MCMC_settings$chain_properties$record_every)

  logpost_chain[1] <- lposterior_total(curr_aug_dat,
                                       curr_theta,
                                       obs_dat,
                                       hyperparameters,
                                       index_dates,
                                       range_dates)

  # Track acceptance counts for all move types
  n_accepted_D_moves <- 0
  n_proposed_D_moves <- 0

  n_accepted_E_moves <- 0
  n_proposed_E_moves <- 0

  n_accepted_swapE_moves <- 0
  n_proposed_swapE_moves <- 0

  n_accepted_mu_moves <- lapply(seq_len(n_groups),
                                function(g) rep(0, ncol(index_dates[[g]])))
  n_proposed_mu_moves <- n_accepted_mu_moves

  n_accepted_CV_moves <- n_accepted_mu_moves
  n_proposed_CV_moves <- n_accepted_mu_moves

  #----------------------------------------------------------------------------
  # Run the MCMC

  print("... Burnin ...")
  for (k in seq_len(MCMC_settings$chain_properties$n_iter - 1)) {

    output_stuff <- (k >= MCMC_settings$chain_properties$burnin) &
                    (k %% MCMC_settings$chain_properties$record_every) == 0

    if (output_stuff) {
      print(sprintf("... %d / %d ...", k,
                    MCMC_settings$chain_properties$n_iter))
    }

    # Move some of the D_i (augmented event dates) ----------------------------
    if (MCMC_settings$moves_switch$D_on) {

      # Loop over each group
      for (g in seq_len(n_groups)) {

        # Loop over each date column in that group
        for(j in seq_len(ncol(curr_aug_dat$D[[g]]))) {

          # propose moves for only a certain fraction of dates
          to_update <- sample(seq_len(nrow(obs_dat[[g]])),
            round(nrow(obs_dat[[g]]) *
                  MCMC_settings$moves_options$fraction_Di_to_update))

          n_groups_to_update <- floor(
            length(to_update) /
            MCMC_settings$moves_options$move_D_by_groups_of_size
            )

          for(i in seq_len(n_groups_to_update)) {

            indices <- to_update[
              MCMC_settings$moves_options$move_D_by_groups_of_size * (i - 1) +
                (seq_len(MCMC_settings$moves_options$move_D_by_groups_of_size))
              ]
            
            tmp <- move_Di(indiv_idx,
                           g,
                           j,
                           curr_aug_dat,
                           curr_theta,
                           obs_dat,
                           hyperparameters,
                           index_dates,
                           range_dates)

            n_proposed_D_moves <- n_proposed_D_moves + 1
            n_accepted_D_moves <- n_accepted_D_moves + tmp$accept

            # if accepted move, update accordingly
            if(tmp$accept == 1) curr_aug_dat <- tmp$new_aug_dat
          }
        }
      }
    }

    # move some of the E_i (error indicators) ---------------------------------
    if (MCMC_settings$moves_switch$E_on) {

      # Loop over each group
      for (g in seq_len(n_groups)) {

        # Loop over each date column in each group
        for(j in seq_len(ncol(curr_aug_dat$E[[g]]))) {

          # proposing moves for only a certain fraction of dates
          to_update <- sample(seq_len(nrow(obs_dat[[g]])),
            round(nrow(obs_dat[[g]]) *
                  MCMC_settings$moves_options$fraction_Ei_to_update))

          n_groups_to_update <- length(to_update)

          for (i in seq_len(n_groups_to_update)) {
            tmp <- move_Ei(to_update[i],
                            g,
                            j,
                            curr_aug_dat,
                            curr_theta,
                            obs_dat,
                            hyperparameters,
                            index_dates,
                            range_dates)

            n_proposed_E_moves <- n_proposed_E_moves + 1
            n_accepted_E_moves <- n_accepted_E_moves + tmp$accept

            # if accepted move, update accordingly
            if (tmp$accept == 1) curr_aug_dat <- tmp$new_aug_dat
          }
        }
      }
    }

    # swap eligible E values (i.e. 0 <-> 1)
    if (MCMC_settings$moves_switch$swapE_on) {

      # Loop over each group
      for (g in seq_len(n_groups)) {

        candidates_for_swap <- find_Eis_to_swap(g, curr_aug_dat)

        for (i in candidates_for_swap) {

          tmp <- swap_Ei(i,
                         g,
                         curr_aug_dat,
                         curr_theta,
                         obs_dat,
                         hyperparameters,
                         index_dates,
                         range_dates)

          n_proposed_swapE_moves <- n_proposed_swapE_moves + 1
          n_accepted_swapE_moves <- n_accepted_swapE_moves + tmp$accept
          # if accepted move, update accordingly
          if (tmp$accept == 1) curr_aug_dat <- tmp$new_aug_dat
        }
      }
    }

    # Update zeta using Gibbs sampling ----------------------------------------
    if (MCMC_settings$moves_switch$zeta_on) {
      tmp <- move_zeta_gibbs(curr_aug_dat, curr_theta, hyperparameters)
      curr_theta <- tmp$new_theta
    }

    # Move mu and CV using log-normal proposals -------------------------------
    for (param in c("mu", "CV")) {
      if (MCMC_settings$moves_switch[[paste0(param,"_on")]]) {
        for (g in seq_len(n_groups)) {
          for (j in seq_len(ncol(index_dates[[g]]))) { #for (j in seq(2, ncol(curr_aug_dat$D[[g]]))) {
            tmp <- move_lognormal(
              what = param,
              group_idx = g,
              delay_idx = j, # delay_idx = j - 1
              sdlog = MCMC_settings$moves_options[[paste0("sdlog_", param)]][[g]][[j]], # [[j - 1]]
              aug_dat = curr_aug_dat,
              curr_theta = curr_theta,
              obs_dat = obs_dat,
              hyperparameters = hyperparameters,
              index_dates = index_dates)
            n_proposed <- get(paste0("n_proposed_", param, "_moves"))
            n_accepted <- get(paste0("n_accepted_", param, "_moves"))
            n_proposed[[g]][j] <- n_proposed[[g]][j] + 1 # n_proposed[[g]][j - 1] <- n_proposed[[g]][j - 1] + 1
            if (tmp$accept == 1) {
              n_accepted[[g]][j] <- n_accepted[[g]][j] + 1 # n_accepted[[g]][j - 1] <- n_accepted[[g]][j - 1] + 1
              curr_theta <- tmp$new_theta
            }
            assign(paste0("n_proposed_", param, "_moves"), n_proposed)
            assign(paste0("n_accepted_", param, "_moves"), n_accepted)
          }
        }
      }
    }

    # Record parameter values and corresponding posterior after all moves -----
    if (output_stuff) {
      idx <- (k - MCMC_settings$chain_properties$burnin) /
              MCMC_settings$chain_properties$record_every + 1
      theta_chain[[idx]] <- curr_theta
      aug_dat_chain[[idx]] <- curr_aug_dat
      logpost_chain[idx] <- lposterior_total(curr_aug_dat,
                                             curr_theta,
                                             obs_dat,
                                             hyperparameters,
                                             index_dates,
                                             range_dates)
      # Previous note: CONSIDER DOING THIS USING SAPPLY AFTER THE WHOLE THING
    }
  }

  #----------------------------------------------------------------------------
  # Compute acceptance probabilities

  accept_prob <- list(
    D_moves = n_accepted_D_moves / n_proposed_D_moves,
    E_moves = n_accepted_E_moves / n_proposed_E_moves,
    mu_moves = lapply(seq_len(n_groups),
                      function(g) {
                        n_accepted_mu_moves[[g]] / n_proposed_mu_moves[[g]]
                        }),
    CV_moves = lapply(seq_len(n_groups),
                      function(g) {
                        n_accepted_CV_moves[[g]] / n_proposed_CV_moves[[g]]
                        }),
    zeta_moves = 1 # always accepted (Gibbs)
    )

  #----------------------------------------------------------------------------
  # Return list of outputs of interest

  res <- list(theta_chain = theta_chain,
              aug_dat_chain = aug_dat_chain,
              logpost_chain = logpost_chain,
              accept_prob = accept_prob,
              index_dates = index_dates)

  return(res)
}

#' Compute correlation between the MCMC chains of mean and CV for each delay
#'
#' @description
#' This function computes Pearson correlation coefficients between the posterior
#'  samples of mu and CV for each delay distribution. Optionally, it produces
#'   scatter plots to visually inspect correlations.
#'
#'
#' @param MCMCres The output of function \code{\link{RunMCMC}}.
#' @param plot A boolean indicating whether to plot the correlations.
#' @param group_labels Optional character vector of group names for more
#'  informative plots. If \code{NULL}, group index will be used.
#' @param date_labels Optional list of character vectors giving the names of
#'  each date column per group. This is used to automatically generate delay
#'   labels. If \code{NULL}, delay index will be used.
#'
#' @return A list of correlation test results (obtained using
#'  \code{\link{cor.test}}), one per group and delay, assessing correlations
#'   between the posterior mean and posterior CV for each delay.
#' @import graphics
#' @importFrom stats cor.test
#' @import ggplot2
#' @import dplyr
#' @export
#' @examples
#' # Simulate data to use
#' n_groups <- 4
#' index_dates <- list(
#'   matrix(c(1, 2), nrow = 2),
#'   cbind(c(1, 2), c(1, 3)),
#'   cbind(c(1, 2), c(2, 3),c(1, 4)),
#'   cbind(c(1, 2), c(2, 3), c(1, 4))
#' )
#' theta <- list(
#'   mu = list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13)),
#'   CV = list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5)),
#'   prop_missing_data = 0.2,
#'   zeta = 0.05
#'  )
#'
#' n_per_group <- rep(10, n_groups)
#' range_dates <- c(0, 30)
#'
#' simul_dat <- simul_true_data(theta, n_per_group, range_dates, index_dates,
#'                              simul_error = TRUE)
#' obs_dat <- simul_dat$obs_dat
#'
#' # Set up hyperparameters
#' hyperparameters <- list(
#'     shape1_prob_error = 3,
#'     shape2_prob_error = 12,
#'     mean_mean_delay = 10,
#'     mean_CV_delay = 10)
#'
#' # Set up MCMC
#' MCMC_settings <- list(
#' moves_switch = list(D_on = TRUE, E_on = TRUE, swapE_on = TRUE,
#'                     mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
#'   moves_options = list(
#'     fraction_Di_to_update = 1 / 10,
#'     move_D_by_groups_of_size = 1,
#'     fraction_Ei_to_update = 1 / 10,
#'     sdlog_mu = list(
#'       0.05,
#'       c(0.15, 0.15),
#'       c(0.15, 0.15, 0.15),
#'       c(0.25, 0.25, 0.25)
#'     ),
#'     sdlog_CV = list(
#'       0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
#'   ),
#'   init_options = list(
#'     mindelay = 0,
#'     maxdelay = 20
#'   ),
#'   chain_properties = list(
#'     n_iter = 500,
#'     burnin = 50,
#'     record_every = 10
#'   )
#' )
#'
#' # Run MCMC
#' MCMC_result <- RunMCMC(obs_dat,
#'                        MCMC_settings,
#'                        hyperparameters,
#'                        index_dates)
#'
#' compute_correlations_mu_CV(MCMCres = MCMC_result,
#'                            plot = TRUE,
#'                            group_labels = c("Community-alive",
#'                                             "Community-dead",
#'                                             "Hospitalised-alive",
#'                                             "Hospitalised-dead"),
#'                           date_labels = list(
#'                               c("Onset", "Report"),
#'                               c("Onset", "Death", "Report"),
#'                               c("Onset", "Hosp", "Disch", "Report"),
#'                               c("Onset", "Hosp", "Death", "Report")
#'                             ))
compute_correlations_mu_CV <- function(MCMCres,
                                       plot = TRUE,
                                       group_labels = NULL,
                                       date_labels = NULL) {

  iterations <- seq_along(MCMCres$theta_chain)
  n_groups <- length(MCMCres$theta_chain[[1]]$mu)
  n_delays <- sapply(MCMCres$theta_chain[[1]]$mu, length)

  # default group labels
  if (is.null(group_labels)) group_labels <- paste("Group", seq_len(n_groups))

  # create delay labels
  if (!is.null(date_labels)) {
    delay_labels <- lapply(seq_along(MCMCres$index_dates), function(g) {
      idx_mat <- MCMCres$index_dates[[g]]
      labels <- date_labels[[g]]
      apply(idx_mat, 2, function(col) paste0(labels[col[1]], " to ", labels[col[2]]))
    })
  } else {
    # default delay labels
    delay_labels <- lapply(n_delays, function(n) paste("Delay", seq_len(n)))
  }

  cor_mu_CV <- list()
  df_list <- list()

  for (g in seq_len(n_groups)) {
    cor_mu_CV[[g]] <- list()
    for (d in seq_len(n_delays[g])) {
      mu_vals <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[g]][d])
      CV_vals <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[g]][d])

      cor_mu_CV[[g]][[d]] <- cor.test(mu_vals, CV_vals)

      df_list[[length(df_list) + 1]] <- data.frame(
        mu = mu_vals,
        CV = CV_vals,
        group = group_labels[g],
        delay = delay_labels[[g]][d]
      )
    }
  }

  # plot
  if (plot) {

    full_df <- bind_rows(df_list) %>%
      mutate(panel_label = paste0(group, ": ", delay))

    cor_plot <- ggplot(full_df, aes(x = mu, y = CV)) +
      geom_point(alpha = 0.5, size = 1.2) +
      facet_wrap(~ panel_label) +
      theme_minimal(base_size = 12) +
      labs(x = "Mu", y = "CV") +
      theme(
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "grey", fill = NA),
        plot.title = element_text(hjust = 0.5)
        )
    print(cor_plot)
  }

  return(cor_mu_CV)
}



#' Compute autocorrelation for each parameter in the MCMC chains
#'
#' @param MCMCres The output of function \code{\link{RunMCMC}}.
#' @param group_labels Optional character vector of group names for more
#'  informative plots. If \code{NULL}, group index will be used.
#' @param date_labels A list of character vectors giving the names of each date
#'  column per group. This is used to automatically generate delay labels. If
#'   \code{NULL}, delay index will be used.
#' @return A list of autocorrelation results and plots for each parameter: mu,
#'  CV and zeta (obtained using the function \code{\link{acf}}).
#' @import graphics
#' @importFrom stats acf
#' @export
#' @examples
#' # Simulate data to use
#' n_groups <- 4
#' index_dates <- list(
#'   matrix(c(1, 2), nrow = 2),
#'   cbind(c(1, 2), c(1, 3)),
#'   cbind(c(1, 2), c(2, 3),c(1, 4)),
#'   cbind(c(1, 2), c(2, 3), c(1, 4))
#' )
#' theta <- list(
#'   mu = list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13)),
#'   CV = list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5)),
#'   prop_missing_data = 0.2,
#'   zeta = 0.05
#'  )
#'
#' n_per_group <- rep(10, n_groups)
#' range_dates <- c(0, 30)
#'
#' simul_dat <- simul_true_data(theta, n_per_group, range_dates, index_dates,
#'                              simul_error = TRUE)
#' obs_dat <- simul_dat$obs_dat
#'
#' # Set up hyperparameters
#' hyperparameters <- list(
#'     shape1_prob_error = 3,
#'     shape2_prob_error = 12,
#'     mean_mean_delay = 10,
#'     mean_CV_delay = 10)
#'
#' # Set up MCMC
#' MCMC_settings <- list(
#' moves_switch = list(D_on = TRUE, E_on = TRUE, swapE_on = TRUE,
#'                     mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
#'   moves_options = list(
#'     fraction_Di_to_update = 1 / 10,
#'     move_D_by_groups_of_size = 1,
#'     fraction_Ei_to_update = 1 / 10,
#'     sdlog_mu = list(
#'       0.05,
#'       c(0.15, 0.15),
#'       c(0.15, 0.15, 0.15),
#'       c(0.25, 0.25, 0.25)
#'     ),
#'     sdlog_CV = list(
#'       0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
#'   ),
#'   init_options = list(
#'     mindelay = 0,
#'     maxdelay = 20
#'   ),
#'   chain_properties = list(
#'     n_iter = 500,
#'     burnin = 50,
#'     record_every = 10
#'   )
#' )
#'
#' # Run MCMC
#' MCMC_result <- RunMCMC(obs_dat,
#'                        MCMC_settings,
#'                        hyperparameters,
#'                        index_dates)
#'
#' compute_autocorr(MCMCres = MCMC_result,
#'                  group_labels = c("Community-alive",
#'                                   "Community-dead",
#'                                   "Hospitalised-alive",
#'                                   "Hospitalised-dead"),
#'                  date_labels = list(
#'                    c("Onset", "Report"),
#'                    c("Onset", "Death", "Report"),
#'                    c("Onset", "Hosp", "Disch", "Report"),
#'                    c("Onset", "Hosp", "Death", "Report")
#'                  ))
compute_autocorr <- function(MCMCres,
                             group_labels = NULL,
                             date_labels = NULL) {

  autocorr <- list(mu = list(), CV = list())

  n_groups <- length(MCMCres$index_dates)
  n_delays <- sapply(MCMCres$index_dates, ncol)
  iterations <- seq_len(length(MCMCres$theta_chain))

  # Create delay labels if date_labels are provided
  if (!is.null(date_labels)) {
    delay_labels <- lapply(seq_along(MCMCres$index_dates), function(g) {
      idx_mat <- MCMCres$index_dates[[g]]
      labels <- date_labels[[g]]
      apply(idx_mat, 2, function(col) paste0(labels[col[1]], " to ", labels[col[2]]))
    })
  } else {
    delay_labels <- lapply(n_delays, function(n) paste0("Delay ", seq_len(n)))
  }

  if (is.null(group_labels)) {
    group_labels <- paste("Group", seq_len(n_groups))
  }

  # Estimate plot layout
  total_plots <- sum(n_delays) * 2 + 1  # mu + CV + zeta
  n_cols <- ceiling(sqrt(total_plots))
  n_rows <- ceiling(total_plots / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

  for (g in seq_len(n_groups)) {
    autocorr$mu[[g]] <- list()
    autocorr$CV[[g]] <- list()

    for (d in seq_len(n_delays[g])) {
      mu_chain <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[g]][d])
      CV_chain <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[g]][d])

      mu_label <- sprintf("Mu: %s, %s", group_labels[g], delay_labels[[g]][d])
      CV_label <- sprintf("CV: %s, %s", group_labels[g], delay_labels[[g]][d])

      autocorr$mu[[g]][[d]] <- acf(mu_chain, main = mu_label)
      autocorr$CV[[g]][[d]] <- acf(CV_chain, main = CV_label)
    }
  }

  # zeta
  zeta_chain <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta)
  autocorr$zeta <- acf(zeta_chain, main = "Zeta")

  return(autocorr)
}

#' Computes posterior estimates of parameters from the MCMC chain
#'
#' @param MCMCres Output from \code{\link{RunMCMC}}.
#' @param central A character string specifying what the central estimate
#'  should be (either \code{"median"} or \code{"mean"} posterior)
#' @param CrI A scalar in [0;1] used to compute the posterior credible
#'  intervals. For 95\% credible intervals, use CrI=0.95.
#' @param theta_true Optional list of true parameter values to compare the
#'  output to. If not \code{NULL}, this should contain:
#' \itemize{
#'  \item{\code{mu}: A list of length
#'   \code{n_groups = length(MCMCres$aug_dat_chain[[1]]$D)}. Each element of
#'    \code{mu} should be a scalar or vector giving the mean delay(s) in that
#'     group.}
#'  \item{\code{CV}: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar or vector giving the coefficient of variation
#'    of the delay(s) in that group.}
#'  \item{\code{zeta}: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' The posterior distributions of parameters are then plotted together with
#'  \code{theta_true}.
#' @param plot A boolean specifying whether to generate boxplots of the
#'  posterior estimates.
#' @param group_labels A character vector of length equal to the number of
#'  groups. Each element gives the label corresponding to a group index e.g.
#'   \code{c("Community-alive", "Hospitalised-dead")}.
#' @param date_labels A list of character vectors giving the names of each date
#'  column per group. This is used to automatically generate delay labels.
#'
#' @return A list with posterior estimates and optional plots:
#' \itemize{
#'  \item{\code{logpost}: A vector of three values, containing the central
#'   log-posterior estimate (mean or median), lower bound of the credible
#'    interval and upper bound of the credible interval.}
#'  \item{\code{theta}: A named list giving posterior parameter estimates for:
#'  \itemize{
#'  \item{\code{mu}: A list of length equal to the number of groups. Each
#'   element is a matrix with 3 rows giving the posterior mean/median of the
#'     delay(s) in that group (1st row = central posterior estimate, 2nd and
#'      3rd rows = lower and upper credible interval bounds).}
#'  \item{\code{CV}: Same structure as \code{mu}, giving posterior estimates
#'   of the coefficient of variation for each delay.}
#'  \item{\code{zeta}: A vector of three values in [0;1] giving the posterior
#'   estimate of the probability that, if a data point is not missing, it is
#'    recorded with error (1st value = central posterior estimate, 2nd and 3rd
#'     values = lower and upper credible interval bounds).}
#'  }}
#' }
#' If \code{plot = TRUE}, a summary plot of the posterior distributions is also
#'  displayed. This includes boxplots for the log-posterior, zeta, and the mu
#'   and cv of delays for each group.
#'
#' @import graphics
#' @import ggplot2
#' @import patchwork
#' @importFrom colorspace scale_fill_discrete_qualitative qualitative_hcl
#' @export
#'
#' @seealso \code{\link{RunMCMC}}
#'
#' @examples
#' # Simulate data to use
#' n_groups <- 4
#' index_dates <- list(
#'   matrix(c(1, 2), nrow = 2),
#'   cbind(c(1, 2), c(1, 3)),
#'   cbind(c(1, 2), c(2, 3),c(1, 4)),
#'   cbind(c(1, 2), c(2, 3), c(1, 4))
#' )
#' theta <- list(
#'   mu = list(5, c(6, 7), c(8, 9, 10), c(11, 12, 13)),
#'   CV = list(0.5, c(0.5, 0.5), c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5)),
#'   prop_missing_data = 0.2,
#'   zeta = 0.05
#'  )
#'
#' n_per_group <- rep(10, n_groups)
#' range_dates <- c(0, 30)
#'
#' simul_dat <- simul_true_data(theta, n_per_group, range_dates, index_dates,
#'                              simul_error = TRUE)
#' obs_dat <- simul_dat$obs_dat
#'
#' # Set up hyperparameters
#' hyperparameters <- list(
#'     shape1_prob_error = 3,
#'     shape2_prob_error = 12,
#'     mean_mean_delay = 10,
#'     mean_CV_delay = 10)
#'
#' # Set up MCMC
#' MCMC_settings <- list(
#' moves_switch = list(D_on = TRUE, E_on = TRUE, swapE_on = TRUE,
#'                     mu_on = TRUE, CV_on = TRUE, zeta_on = TRUE),
#'   moves_options = list(
#'     fraction_Di_to_update = 1 / 10,
#'     move_D_by_groups_of_size = 1,
#'     fraction_Ei_to_update = 1 / 10,
#'     sdlog_mu = list(
#'       0.05,
#'       c(0.15, 0.15),
#'       c(0.15, 0.15, 0.15),
#'       c(0.25, 0.25, 0.25)
#'     ),
#'     sdlog_CV = list(
#'       0.25, c(0.25, 0.25), c(0.25, 0.25, 0.25), c(0.25, 0.25, 0.25))
#'   ),
#'   init_options = list(
#'     mindelay = 0,
#'     maxdelay = 20
#'   ),
#'   chain_properties = list(
#'     n_iter = 500,
#'     burnin = 50,
#'     record_every = 10
#'   )
#' )
#'
#' # Run MCMC
#' MCMC_result <- RunMCMC(obs_dat,
#'                        MCMC_settings,
#'                        hyperparameters,
#'                        index_dates)
#'
#' # Get summary
#' get_param_posterior_estimates(MCMCres = MCMC_result,
#'                               central = "mean",
#'                               CrI = 0.95,
#'                               theta_true = NULL,
#'                               plot = TRUE,
#'                               group_labels = c("Community-alive",
#'                                                "Community-dead",
#'                                                "Hospitalised-alive",
#'                                                "Hospitalised-dead"),
#'                               date_labels = list(
#'                                 c("Onset", "Report"),
#'                                 c("Onset", "Death", "Report"),
#'                                 c("Onset", "Hosp", "Disch", "Report"),
#'                                 c("Onset", "Hosp", "Death", "Report")
#'                               ))
#'
#'
get_param_posterior_estimates <- function(MCMCres,
                                          central = c("median", "mean"),
                                          CrI = 0.95,
                                          theta_true = NULL,
                                          plot = TRUE,
                                          group_labels = NULL,
                                          date_labels = NULL) {

  # checks -------------------------------------------------------------------

  if (is.null(group_labels)) stop("Supply group_labels (character vector)")
  if (length(group_labels) != length(MCMCres$aug_dat_chain[[1]]$D)) {
    stop("group_labels length must match the number of groups")
  }
  if (is.null(date_labels)) stop("Supply date_labels (list of character vectors)")
  if (is.null(MCMCres$index_dates)) stop("MCMCres must contain index_dates")

  # create delay labels using index_dates and date_labels --------------------
  generate_delay_labels <- function(index_dates, date_labels) {
    delay_labels <- list()
    delay_orders <- list()

    for (g in seq_along(index_dates)) {
      idx_mat <- index_dates[[g]]
      labels <- date_labels[[g]]

      delays <- apply(idx_mat, 2, function(col) paste0(labels[col[1]], " to ", labels[col[2]]))

      # Ordering by end date (col[2]) then start date (col[1])
      ordering <- order(idx_mat[2, ], idx_mat[1, ])

      delay_labels[[g]] <- delays[ordering]
      delay_orders[[g]] <- ordering
    }

    list(labels = delay_labels, order = delay_orders)
  }

  index_dates <- MCMCres$index_dates
  delay_info <- generate_delay_labels(index_dates, date_labels)
  delay_labels <- delay_info$labels
  delay_order <- delay_info$order

  delay_levels_df <- do.call(rbind, lapply(seq_along(group_labels), function(g) {
    data.frame(
      group = group_labels[g],
      delay = delay_labels[[g]],
      delay_order = delay_order[[g]]
    )
  }))

  iterations <- seq_len(length(MCMCres$theta_chain))
  output <- list()

  # Log posterior summary ----------------------------------------------------
  output$logpost <- c(
    get(central)(MCMCres$logpost_chain),
    quantile(MCMCres$logpost_chain, c((1 - CrI) / 2, CrI + (1 - CrI) / 2))
  )

  # Get parameter estimates --------------------------------------------------
  output$theta <- list(mu = list(), CV = list())
  if (plot) plot_data <- list()

  for (group_idx in seq_along(group_labels)) {
    for (param in c("mu", "CV")) {
      values <- lapply(seq_len(ncol(index_dates[[group_idx]])), function(j) {
        sapply(iterations, function(e) MCMCres$theta_chain[[e]][[param]][[group_idx]][j])
      })

      est_matrix <- sapply(seq_along(values), function(j) {
        c(get(central)(values[[j]]),
          quantile(values[[j]], c((1 - CrI) / 2, 1 - (1 - CrI) / 2)))
      })

      output$theta[[param]][[group_idx]] <- est_matrix

      if (plot) {
        plot_data[[length(plot_data) + 1]] <- do.call(
          rbind, lapply(seq_along(values), function(j) {
          data.frame(
            value = values[[j]],
            delay = delay_labels[[group_idx]][j],
            group = group_labels[group_idx],
            param = param
          )
        })
        )
      }
    }
  }

  # zeta (scalar, not per group) ---------------------------------------------
  zeta_chain <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta)
  output$theta$zeta <- c(get(central)(zeta_chain),
                         quantile(zeta_chain, c((1 - CrI) / 2, 1 - (1 - CrI) / 2)))

  if (plot) {
    central_label <- central
    ylab_central <- ifelse(central_label == "median",
                           "Posterior Median",
                           "Posterior Mean")

    # plot log posterior -----------------------------------------------------
    df_logpost <- data.frame(
      value = MCMCres$logpost_chain,
      param = "Log Posterior"
    )

    p1 <- ggplot(df_logpost, aes(x = param, y = value)) +
      geom_boxplot(width = 0.5, fill = "gray90") +
      labs(x = NULL, y = ylab_central) +
      facet_wrap(~ param) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = "grey"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(colour = "grey", fill = NA),
        axis.title.y = element_text(margin = margin(r = 10))
      )

    # plot zeta --------------------------------------------------------------
    df_zeta <- data.frame(value = zeta_chain, param = "zeta")

    p2 <- ggplot(df_zeta, aes(x = param, y = value)) +
      geom_boxplot(width = 0.5, fill = "gray90") +
      labs(x = NULL, y = ylab_central) +
      facet_wrap(~ param) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = "grey"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(colour = "grey", fill = NA),
        axis.title.y = element_text(margin = margin(r = 10))
      )

    # Add reference line for true zeta if supplied
    if (!is.null(theta_true) && !is.null(theta_true$zeta)) {
      p2 <- p2 + geom_hline(yintercept = theta_true$zeta,
                            linetype = "dashed", color = "black")
    }

    # plot mu delays ---------------------------------------------------------
    df_all <- do.call(rbind, plot_data)
    df_params <- df_all[df_all$param %in% c("mu", "CV"), ]

    df_mu <- subset(df_params, param == "mu")

    df_mu <- df_mu %>%
      left_join(delay_levels_df, by = c("group", "delay")) %>%
      mutate(delay = factor(delay, levels = unique(delay[order(delay_order)]))) %>%
      select(-delay_order)

    # Consistent colours
    all_delays <- sort(unique(df_mu$delay))
    delay_colours <- qualitative_hcl(length(all_delays), palette = "Dynamic")
    names(delay_colours) <- all_delays

    p3 <- ggplot(df_mu, aes(x = delay, y = value, fill = delay)) +
      geom_boxplot(width = 0.5, colour = "black", alpha = 0.8) +
      facet_wrap(~ group, scales = "free_x", nrow = 1) +
      labs(x = NULL, y = paste(ylab_central, "of Mu delays")) +
      scale_fill_manual(values = delay_colours) +
      scale_colour_manual(values = delay_colours) +
      theme_minimal(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "grey90", colour = "grey"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(colour = "grey", fill = NA),
        axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "none"
      )

    # Add reference lines for true mu values if supplied
    if (!is.null(theta_true) && !is.null(theta_true$mu)) {
      df_mu_true <- do.call(rbind, lapply(seq_along(theta_true$mu), function(g) {
        data.frame(
          delay = delay_labels[[g]],
          group = group_labels[g],
          true_value = theta_true$mu[[g]]
        )
      }))
      p3 <- p3 + geom_hline(
        data = df_mu_true,
        aes(yintercept = true_value, group = delay, colour = delay),
        linetype = "dashed"
      )
    }

    # plot cv delays ---------------------------------------------------------
    df_cv <- subset(df_params, param == "CV")

    df_cv <- df_cv %>%
      left_join(delay_levels_df, by = c("group", "delay")) %>%
      mutate(delay = factor(delay, levels = unique(delay[order(delay_order)]))) %>%
      select(-delay_order)

    p4 <- ggplot(df_cv, aes(x = delay, y = value, fill = delay)) +
      geom_boxplot(width = 0.5, colour = "black", alpha = 0.8) +
      facet_wrap(~ group, scales = "free_x", nrow = 1) +
      labs(x = NULL, y = paste(ylab_central, "of CV delays")) +
      scale_fill_manual(values = delay_colours) +
      scale_colour_manual(values = delay_colours) +
      scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
      theme_minimal(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "grey90", colour = "grey"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(colour = "grey", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "none"
      )

    # Add reference lines for true CV values if supplied
    if (!is.null(theta_true) && !is.null(theta_true$CV)) {
      df_cv_true <- do.call(rbind, lapply(seq_along(theta_true$CV), function(g) {
        data.frame(
          delay = delay_labels[[g]],
          group = group_labels[g],
          true_value = theta_true$CV[[g]]
        )
      }))
      p4 <- p4 + geom_hline(
        data = df_cv_true,
        aes(yintercept = true_value, group = delay, colour = delay),
        linetype = "dashed"
      )
    }

    # plot together ----------------------------------------------------------
    left_column <- p1 / p2
    right_column <- p3 / p4
    combined_plot <- (left_column | plot_spacer() | right_column) +
      plot_layout(widths = c(2, 0.2, 8))

    print(combined_plot)
  }
  return(output)
}
