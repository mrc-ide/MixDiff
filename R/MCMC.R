#' Runs the MCMC estimation procedure.
#' 
#' @param obs_dat A list of observed data, in the format of the first element
#'  (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param MCMC_settings A list of settings to be used for running the MCMC, see
#'  details.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider
#'  in the estimation, see details.
#' @details \code{MCMC_settings} should be a list containing:
#' \itemize{
#'  \item{\code{moves_switch}}{: A list of booleans (D_on ,E_on, mu_on, CV_on,
#'   zeta_on) stating whether each parameter/augmented data should be moved in
#'    the procedure or not.}
#'  \item{\code{moves_options}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{fraction_Di_to_update}}{: The fraction of augmented dates to be
#'   updated at each iteration of the MCMC.}
#'  \item{\code{move_D_by_groups_of_size}}{: The number of augmented dates to be
#'   updated simultaneously in each group.}
#'  \item{\code{fraction_Ei_to_update}}{: The fraction of indicators of whether
#'   observed dates are erroneous to be updated at each iteration of the MCMC.}
#'  \item{\code{sdlog_mu}}{: The standard deviations to be used for proposing
#'   moves of the mean delays. This should be a list of length
#'    \code{n_groups = length(obs_dat)}. Each element in the list should be a
#'     vector with length given by the numbers of delays to be considered in
#'      this group.}
#'  \item{\code{sdlog_CV}}{: The standard deviations to be used for proposing
#'   moves of the CV of delays. This should be a list of length
#'    \code{n_groups = length(obs_dat)}. Each element in the list should be a
#'     vector with length given by the numbers of delays to be considered in
#'      this group.}
#'  }
#'  }
#'  \item{\code{init_options}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{mindelay}}{: The minimum delay, below which dates are considered
#'   incompatile with one another at the initialisation stage of the MCMC.}
#'  \item{\code{maxdelay}}{: The maximum delay, above which dates are considered
#'   incompatile with one another at the initialisation stage of the MCMC.}
#'  \item{\code{record_every}}{: A number indicating, after the burnin, every
#'   how many iterations outputs should be recorded.}
#'  }
#'  }
#'  \item{\code{chain_properties}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{n_iter}}{: The total number of iteration of MCMC to run.}
#'  \item{\code{burnin}}{: The number of initial iterations to consider as the
#'   burnin period - no output is recorded for these initial MCMC iterations.}
#'  \item{\code{record_every}}{: A number indicating, after the burnin, every
#'   how many iterations outputs should be recorded.}
#'  }
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
#'  
#'  theta_chain = theta_chain, aug_dat_chain = aug_dat_chain, logpost_chain =
#'   logpost_chain, accept_prob = accept_prob
#'  \item{\code{theta_chain}}{: a list of parameters, at each recorded step of
#'   the MCMC chain}
#'  \item{\code{aug_dat_chain}}{: a list of augmented data, at each recorded
#'   step of the MCMC chain}
#'  \item{\code{logpost_chain}}{: a vector of values of the log posterior at
#'   each recorded step of the MCMC chain}
#'  \item{\code{accept_prob}}{: A list of the proababilities of acceptance for
#'   each parameter across all MCMC iterations}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
RunMCMC <- function(obs_dat, 
                    MCMC_settings,
                    hyperparameters,
                    index_dates) {
  
  n_dates <- sapply(obs_dat, ncol )
  n_groups <- length(n_dates)
  
  ###############################################
  ### define augmented data to be used for initialisation of the chain ###
  ###############################################
  
  check_MCMC_settings(MCMC_settings, index_dates)
  
  ###############################################
  ### define augmented data to be used for initialisation of the chain ###
  ###############################################
  
  aug_dat <- initialise_aug_data(obs_dat,
                                 compute_index_dates_order(index_dates),
                                 MCMC_settings)
  
  ###############################################
  ### define parameters to be used for initialisation of the chain ###
  ###############################################
  
  theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)
  
  ###############################################
  ### Initalise the MCMC chains ###
  ###############################################
  
  range_dates <- find_range(obs_dat)
  
  # to store param values
  curr_theta <- theta
  theta_chain <- list()
  theta_chain[[1]] <- curr_theta
  
  # to store augmented data values
  curr_aug_dat <- aug_dat
  aug_dat_chain <- list()
  aug_dat_chain[[1]] <- curr_aug_dat
  
  logpost_chain <- rep(NA,
                       (MCMC_settings$chain_properties$n_iter -
                          MCMC_settings$chain_properties$burnin) /
                         MCMC_settings$chain_properties$record_every)
  
  logpost_chain[1] <- lposterior_total(curr_aug_dat,
                                       curr_theta,
                                       obs_dat,
                                       hyperparameters,
                                       index_dates,
                                       range_dates)
  
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
  
  ###############################################
  ### Run the MCMC ###
  ###############################################
  
  print("... Burnin ...")
  
  for (k in seq_len(MCMC_settings$chain_properties$n_iter - 1)) {
    
    output_stuff <- (k >= MCMC_settings$chain_properties$burnin) &
      (k %% MCMC_settings$chain_properties$record_every) == 0
    
    if (output_stuff) {
      print(sprintf("... %d / %d ...", k, MCMC_settings$chain_properties$n_iter))
    }
    
    # move some of the D_i
    if (MCMC_settings$moves_switch$D_on) {
      
      # Loop over each group
      for (g in seq_len(n_groups)) {
        
        # Loop over each date column in that group
        for(j in seq_len(ncol(curr_aug_dat$D[[g]]))) {
          
          # propose moves for only a certain fraction of dates
          to_update <- sample(
            seq_len(nrow(obs_dat[[g]])),
            round(
              nrow(obs_dat[[g]]) * MCMC_settings$moves_options$fraction_Di_to_update
              )
            ) 
          
          n_groups_to_update <- floor(
            length(to_update) / MCMC_settings$moves_options$move_D_by_groups_of_size
            )
          
          for(i in seq_len(n_groups_to_update)) {

            tmp <- move_Di(to_update[MCMC_settings$moves_options$move_D_by_groups_of_size * (i - 1) + (seq_len(MCMC_settings$moves_options$move_D_by_groups_of_size))],
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
            if(tmp$accept == 1) {
              # if accepted move, update accordingly
              curr_aug_dat <- tmp$new_aug_dat
              
              # if accepted move, update zeta
              # tmp <- move_zeta_gibbs(curr_aug_dat,
              #                        curr_theta, 
              #                        obs_dat, 
              #                        hyperparameters) 
              # curr_theta <- tmp$new_theta # always update with new theta (Gibbs sampler)
            }
          }
        }
      }
    }
    
    # move some of the E_i
    if (MCMC_settings$moves_switch$E_on) {
      
      # Loop over each group
      for (g in seq_len(n_groups)) {

        # Loop over each date column in each group
        for(j in seq_len(ncol(curr_aug_dat$E[[g]]))) {
          
          # proposing moves for only a certain fraction of dates
          to_update <- sample(
            seq_len(nrow(obs_dat[[g]])),
            round(nrow(obs_dat[[g]]) * MCMC_settings$moves_options$fraction_Ei_to_update)
            )
          
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
            if (tmp$accept == 1) {
              # if accepted move, update accordingly
              curr_aug_dat <- tmp$new_aug_dat
              
              # if accepted move, update zeta
              # tmp <- move_zeta_gibbs(curr_aug_dat,
              #                        curr_theta,
              #                        obs_dat,
              #                        hyperparameters)
              # curr_theta <- tmp$new_theta # always update with new theta (Gibbs sampler)
            }
          }
        }
      }
    }
    
    # swap the E_is that can be swapped (i.e. where exactly one is =1 and exactly one is =0)
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
          if (tmp$accept == 1) {
            # if accepted move, update accordingly
            curr_aug_dat <- tmp$new_aug_dat
            
            # if accepted move, update zeta
            # tmp <- move_zeta_gibbs(curr_aug_dat,
            #                        curr_theta,
            #                        obs_dat,
            #                        hyperparameters)
            # curr_theta <- tmp$new_theta # always update with new theta (Gibbs sampler)
          }
        }
      }
    }
    
    # move zeta using Gibbs sampler
    #print("Move zeta")
    if (MCMC_settings$moves_switch$zeta_on) {
      tmp <- move_zeta_gibbs(curr_aug_dat,
                             curr_theta, 
                             obs_dat, 
                             hyperparameters)
      # always update with new theta (Gibbs sampler)
      curr_theta <- tmp$new_theta
    }
    
    # move mu
    if (MCMC_settings$moves_switch$mu_on) {
      
      for (g in seq_len(n_groups)) {

        for (j in seq(2, ncol(curr_aug_dat$D[[g]]),1)) {

          tmp <- move_lognormal(what = "mu",
                                g,
                                j - 1,
                                MCMC_settings$moves_options$sdlog_mu[[g]][[j - 1]], 
                                curr_aug_dat,
                                curr_theta, 
                                obs_dat, 
                                hyperparameters,
                                index_dates)
          n_proposed_mu_moves[[g]][j - 1] <- n_proposed_mu_moves[[g]][j - 1] + 1
          n_accepted_mu_moves[[g]][j - 1] <- n_accepted_mu_moves[[g]][j - 1] + tmp$accept
          # if accepted move, update accordingly
          if (tmp$accept == 1) curr_theta <- tmp$new_theta
        }
      }
    }
    
    # move CV
    if (MCMC_settings$moves_switch$CV_on) {
      
      for (g in seq_len(n_groups)) {

        for(j in seq(2, ncol(curr_aug_dat$D[[g]]), 1)) {

          tmp <- move_lognormal(what = "CV",
                                g,
                                j - 1,
                                MCMC_settings$moves_options$sdlog_CV[[g]][[j - 1]], 
                                curr_aug_dat,
                                curr_theta, 
                                obs_dat, 
                                hyperparameters,
                                index_dates)
          n_proposed_CV_moves[[g]][j - 1] <- n_proposed_CV_moves[[g]][j - 1] + 1
          n_accepted_CV_moves[[g]][j - 1] <- n_accepted_CV_moves[[g]][j - 1] + tmp$accept
          # if accepted move, update accordingly
          if (tmp$accept == 1) curr_theta <- tmp$new_theta
        }
      }
    }
    
    # recording value of parameters and corresponding posterior after all moves 
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
                                             range_dates) #### CONSIDER DOING THIS USING SAPPLY AFTER THE WHOLE THING
    }
  }
  
  ###############################################
  ### Compute acceptance probabilities ###
  ###############################################
  
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
    zeta_moves = 1)
  
  ###############################################
  ### Return list of outputs of interest ###
  ###############################################
  
  res <- list(theta_chain = theta_chain,
              aug_dat_chain = aug_dat_chain,
              logpost_chain = logpost_chain,
              accept_prob = accept_prob)
  
  return(res)
}

