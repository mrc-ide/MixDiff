#' Runs the MCMC estimation procedure.
#' 
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @param MCMC_settings A list of settings to be used for running the MCMC, see details.
#' @param hyperpriors A list of hyperpriors: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @param index_dates_order A list containing indications on ordering of dates, see details. #### CONSIDER CALCULATING THIS AUTOMATICALLY FROM index_dates
#' @details \code{MCMC_settings} should be a list containing:
#' \itemize{
#'  \item{\code{moves_switch}}{: A list of booleans (D_on , mu_on, CV_on, zeta_on) stating whether each parameter/augmented data should be moved in the procedure or not.}
#'  \item{\code{moves_options}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{fraction_Di_to_update}}{: The fraction of augmented dates to be updated at each iteration of the MCMC.}
#'  \item{\code{move_D_by_groups_of_size}}{: The number of augmented dates to be updates simultaneously in each group.}
#'  \item{\code{sdlog_mu}}{: The standard deviation to be used for proposing moves of the mean delays.}
#'  \item{\code{sdlog_CV}}{: The standard deviation to be used for proposing moves of the CV of delays.}
#'  }
#'  }
#'  \item{\code{chain_properties}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{n_iter}}{: The total number of iteration of MCMC to run.}
#'  \item{\code{burnin}}{: The number of initial iterations to consider as the burnin period - no output is recorded for these initial MCMC iterations. }
#'  \item{\code{record_every}}{: A number indicating, after the burnin, every how many iterations outputs should be recorded.}
#'  }
#'  }
#' }
#' \code{hyperpriors} should be a list containing:
#' \itemize{
#'  \item{\code{shape1_prob_error}}{: A scalar giving the first shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{shape2_prob_error}}{: A scalar giving the second shape parameter for the beta prior used for parameter \code{theta$zeta}}
#'  \item{\code{mean_mean_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$mu}}
#'  \item{\code{mean_CV_delay}}{: A scalar giving the mean of the exponential prior used for parameter \code{theta$CV}}
#' }
#' \code{index_dates} should be a list of length \code{n_groups=length(obs_dat)}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' The number of columns of index_dates[[k]] should match the length of theta$mu[[k]] and theta$CV[[k]] 
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and c(1, 3), this indicates that theta$mu[[k]] and theta$CV[[k]] are respectively the mean and coefficient of variation of two delays: the first delay being between date 1 and date 2, and the second being between date 1 and date 3. 
#' 
#' \code{index_dates_order} should be a list of length \code{n_groups=length(obs_dat)}. Each element of \code{index_dates_order} should be a matrix with 2 rows and a number of columns corresponding to the delays with order rules for that group. 
#' For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date.
#' Each column specifies a rule saying that the origin date must be before the destination date.  
#' @return A list of two elements:
#'  \itemize{
#'  \item{\code{new_aug_dat}}{: Same as \code{curr_aug_dat} but where the relevant dates have been updated}
#'  \item{\code{accept}}{: A scalar with value 1 if the move was accepted and 0 otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
RunMCMC <- function(obs_dat, 
                    MCMC_settings,
                    hyperpriors,
                    index_dates,
                    index_dates_order) ### CHANGE THIS SO index_dates_order is computed automatically from index_dates
{
  
  n_dates <- sapply(obs_dat, ncol )
  n_groups <- length(n_dates)
  
  ###############################################
  ### define augmented data to be used for initialisation of the chain ###
  ###############################################
  
  aug_dat <- initialise_aug_data(obs_dat, index_dates_order)
  
  ###############################################
  ### define parameters to be used for initialisation of the chain ###
  ###############################################
  
  index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
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
  
  logpost_chain <- rep(NA, (MCMC_settings$chain_properties$n_iter - MCMC_settings$chain_properties$burnin) / MCMC_settings$chain_properties$record_every)
  logpost_chain[1] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, hyperpriors, index_dates, range_dates)
  
  n_accepted_D_moves <- 0
  n_proposed_D_moves <- 0
  
  n_accepted_mu_moves <- 0
  n_proposed_mu_moves <- 0
  
  n_accepted_CV_moves <- 0
  n_proposed_CV_moves <- 0
  
  ###############################################
  ### Run the MCMC ###
  ###############################################
  
  print("... Burnin ...")
  
    for(k in seq_len(MCMC_settings$chain_properties$n_iter-1))
    {
      output_stuff <- (k>=MCMC_settings$chain_properties$burnin) & (k %% MCMC_settings$chain_properties$record_every)==0
      
      if(output_stuff)
      {
        print(sprintf("... %d / %d ...", k, MCMC_settings$chain_properties$n_iter))
      }
      
      # move some of the D_i
      if(MCMC_settings$moves_switch$D_on)
      {
        for(g in seq_len(n_groups))
        {
          for(j in seq_len(ncol(curr_aug_dat$D[[g]])))
          {
            to_update <- sample(1:nrow(obs_dat[[g]]), round(nrow(obs_dat[[g]])*MCMC_settings$moves_options$fraction_Di_to_update)) # proposing moves for only a certain fraction of dates
            n_10_to_update <- floor(length(to_update) / MCMC_settings$moves_options$move_D_by_groups_of_size)
            for(i in seq_len(length(n_10_to_update)))
            {
              tmp <- move_Di (to_update[MCMC_settings$moves_options$move_D_by_groups_of_size*(i-1)+(seq_len(MCMC_settings$moves_options$move_D_by_groups_of_size))], g, j, 
                              curr_aug_dat,
                              curr_theta, 
                              obs_dat, 
                              hyperpriors, 
                              index_dates,
                              range_dates) 
              n_proposed_D_moves <- n_proposed_D_moves + 1
              n_accepted_D_moves <- n_accepted_D_moves + tmp$accept
              if(tmp$accept==1) curr_aug_dat <- tmp$new_aug_dat # if accepted move, update accordingly
            }
          }
        }
      }
      
      # move zeta using Gibbs sampler
      if(MCMC_settings$moves_switch$zeta_on)
      {
        tmp <- move_zeta_gibbs(curr_aug_dat,
                               curr_theta, 
                               obs_dat, 
                               hyperpriors) 
        curr_theta <- tmp$new_theta # always update with new theta (Gibbs sampler)
      }
      
     # move mu
      if(MCMC_settings$moves_switch$mu_on)
      {
        for(g in seq_len(n_groups))
        {
          for(j in seq(2,ncol(curr_aug_dat$D[[g]]),1))
          {
            tmp <- move_lognormal(what="mu", g, j-1, MCMC_settings$moves_options$sdlog_mu, 
                                  curr_aug_dat,
                                  curr_theta, 
                                  obs_dat, 
                                  hyperpriors,
                                  index_dates)
            n_proposed_mu_moves <- n_proposed_mu_moves + 1
            n_accepted_mu_moves <- n_accepted_mu_moves + tmp$accept
            if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
          }
        }
      }
      
      # move CV
      if(MCMC_settings$moves_switch$CV_on)
      {
        for(g in seq_len(n_groups))
        {
          for(j in seq(2,ncol(curr_aug_dat$D[[g]]),1))
          {
            tmp <- move_lognormal(what="CV", g, j-1, MCMC_settings$moves_options$sdlog_CV, 
                                  curr_aug_dat,
                                  curr_theta, 
                                  obs_dat, 
                                  hyperpriors,
                                  index_dates)
            n_proposed_CV_moves <- n_proposed_CV_moves + 1
            n_accepted_CV_moves <- n_accepted_CV_moves + tmp$accept
            if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
          }
        }
      }
      
      # recording value of parameters and corresponding posterior after all moves 
      if( output_stuff )
      {
        idx <- (k - MCMC_settings$chain_properties$burnin) / MCMC_settings$chain_properties$record_every+1
        theta_chain[[idx]] <- curr_theta
        aug_dat_chain[[idx]] <- curr_aug_dat
        logpost_chain[idx] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, hyperpriors, index_dates, range_dates) #### CONSIDER DOING THIS USING SAPPLY AFTER THE WHOLE THING
      }
    }
  
  ###############################################
  ### Compute acceptance probabilities ###
  ###############################################
  
  accept_prob <- list(
    D_moves=n_accepted_D_moves / n_proposed_D_moves,
    mu_moves=n_accepted_mu_moves / n_proposed_mu_moves,
    CV_moves=n_accepted_CV_moves / n_proposed_CV_moves,
    zeta_moves=1)
  
  ###############################################
  ### Return lisst of outputs of interest ###
  ###############################################
  
  res <- list(theta_chain=theta_chain, aug_dat_chain=aug_dat_chain, logpost_chain=logpost_chain, accept_prob=accept_prob)
  
  return(res)
}