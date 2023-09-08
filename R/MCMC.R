#' Runs the MCMC estimation procedure.
#'
#' @param obs_dat A list of observed data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}.
#' @param MCMC_settings A list of settings to be used for running the MCMC, see details.
#' @param hyperparameters A list of hyperparameters: see details.
#' @param index_dates A list containing indications on which delays to consider in the estimation, see details.
#' @param verbose logical. If TRUE messages will be printed. Useful for debugging. Defaults to FALSE
#' @details \code{MCMC_settings} should be a list containing:
#' \itemize{
#'  \item{\code{moves_switch}}{: A list of booleans (D_on ,E_on, mu_on, CV_on, zeta_on) stating whether each parameter/augmented data should be moved in the procedure or not.}
#'  \item{\code{moves_options}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{fraction_Di_to_update}}{: The fraction of augmented dates to be updated at each iteration of the MCMC.}
#'  \item{\code{move_D_by_groups_of_size}}{: The number of augmented dates to be updated simultaneously in each group.}
#'  \item{\code{fraction_Ei_to_update}}{: The fraction of indicators of whether observed dates are erroneous to be updated at each iteration of the MCMC.}
#'  \item{\code{sdlog_mu}}{: The standard deviations to be used for proposing moves of the mean delays. This should be a list of length \code{n_groups=length(obs_dat)}. Each element in the list should be a vector with length given by the numbers of delays to be considered in this group.}
#'  \item{\code{sdlog_CV}}{: The standard deviations to be used for proposing moves of the CV of delays. This should be a list of length \code{n_groups=length(obs_dat)}. Each element in the list should be a vector with length given by the numbers of delays to be considered in this group.}
#'  }
#'  }
#'  \item{\code{init_options}}{: A list of the following elements:
#'  \itemize{
#'  \item{\code{mindelay}}{: The minimum delay, below which dates are considered incompatile with one another at the initialisation stage of the MCMC.}
#'  \item{\code{maxdelay}}{: The maximum delay, above which dates are considered incompatile with one another at the initialisation stage of the MCMC.  }
#'  \item{\code{record_every}}{: A number indicating, after the burnin, every how many iterations outputs should be recorded.}
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
#' \code{hyperparameters} should be a list containing:
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
#' @return A list of the following elements:
#'  \itemize{
#'
#'  theta_chain=theta_chain, aug_dat_chain=aug_dat_chain, logpost_chain=logpost_chain, accept_prob=accept_prob
#'  \item{\code{theta_chain}}{: a list of parameters, at each recorded step of the MCMC chain}
#'  \item{\code{aug_dat_chain}}{: a list of augmented data, at each recorded step of the MCMC chain}
#'  \item{\code{logpost_chain}}{: a vector of values of the log posterior at each recorded step of the MCMC chain}
#'  \item{\code{accept_prob}}{: A list of the proababilities of acceptance for each parameter across all MCMC iterations}
#' \item{\code{convergence}}{: a list of logicals. TRUE if the chain has converged i.e., if Gelman-Rubin diagnostic is less than 1.1, FALSE otherwise}
#' }
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
RunMCMC <- function(obs_dat,
                    MCMC_settings,
                    hyperparameters,
                    index_dates, verbose = FALSE)
{

  n_dates <- sapply(obs_dat, ncol )
  n_groups <- length(n_dates)

  ###############################################
  ### define augmented data to be used for initialisation of the chain ###
  ###############################################

  check_MCMC_settings(MCMC_settings, index_dates)

  ###############################################
  ### define augmented data to be used for initialisation of the chain ###
  ###############################################

  aug_dat <- initialise_aug_data(obs_dat, compute_index_dates_order(index_dates), MCMC_settings)

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

  logpost_chain <- rep(NA, (MCMC_settings$chain_properties$n_iter - MCMC_settings$chain_properties$burnin) / MCMC_settings$chain_properties$record_every)
  logpost_chain[1] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, hyperparameters, index_dates, range_dates)

  n_accepted_D_moves <- 0
  n_proposed_D_moves <- 0

  n_accepted_E_moves <- 0
  n_proposed_E_moves <- 0

  n_accepted_swapE_moves <- 0
  n_proposed_swapE_moves <- 0

  n_accepted_mu_moves <- lapply(seq_len(n_groups), function(g) rep(0, ncol(index_dates[[g]]) ))
  n_proposed_mu_moves <- n_accepted_mu_moves

  n_accepted_CV_moves <- n_accepted_mu_moves
  n_proposed_CV_moves <- n_accepted_mu_moves

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
    #print("Move some D_i")
    if(MCMC_settings$moves_switch$D_on)
    {
      for(g in seq_len(n_groups))
      {
        #print(g)
        for(j in seq_len(ncol(curr_aug_dat$D[[g]])))
        {
          #print(j)
          to_update <- sample(seq_len(nrow(obs_dat[[g]])), round(nrow(obs_dat[[g]])*MCMC_settings$moves_options$fraction_Di_to_update)) # proposing moves for only a certain fraction of dates
          n_groups_to_update <- floor(length(to_update) / MCMC_settings$moves_options$move_D_by_groups_of_size)
          for(i in seq_len(n_groups_to_update))
          {
            #print(i)
            tmp <- move_Di (to_update[MCMC_settings$moves_options$move_D_by_groups_of_size*(i-1)+(seq_len(MCMC_settings$moves_options$move_D_by_groups_of_size))], g, j,
                            curr_aug_dat,
                            curr_theta,
                            obs_dat,
                            hyperparameters,
                            index_dates,
                            range_dates)
            n_proposed_D_moves <- n_proposed_D_moves + 1
            n_accepted_D_moves <- n_accepted_D_moves + tmp$accept
            if(tmp$accept==1)
            {
              curr_aug_dat <- tmp$new_aug_dat # if accepted move, update accordingly

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
    #print("Move some E_i")
    if(MCMC_settings$moves_switch$E_on)
    {
      for(g in seq_len(n_groups))
      {
        #print(g)
        for(j in seq_len(ncol(curr_aug_dat$E[[g]])))
        {
          #print(j)
          to_update <- sample(seq_len(nrow(obs_dat[[g]])), round(nrow(obs_dat[[g]])*MCMC_settings$moves_options$fraction_Ei_to_update)) # proposing moves for only a certain fraction of dates
          n_groups_to_update <- length(to_update)
          for(i in seq_len(n_groups_to_update))
          {
            tmp <- move_Ei (to_update[i], g, j,
                            curr_aug_dat,
                            curr_theta,
                            obs_dat,
                            hyperparameters,
                            index_dates,
                            range_dates)
            n_proposed_E_moves <- n_proposed_E_moves + 1
            n_accepted_E_moves <- n_accepted_E_moves + tmp$accept
            if(tmp$accept==1)
            {
              curr_aug_dat <- tmp$new_aug_dat # if accepted move, update accordingly

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
    #print("Swap some E_i")
    if(MCMC_settings$moves_switch$swapE_on)
    {
      for(g in seq_len(n_groups))
      {
        #print(g)
        candidates_for_swap <- find_Eis_to_swap(g, curr_aug_dat)
        #print(candidates_for_swap)
        for(i in candidates_for_swap)
        {
          #print(i)
          tmp <- swap_Ei(i, g,
                  curr_aug_dat,
                  curr_theta,
                  obs_dat,
                  hyperparameters,
                  index_dates,
                  range_dates)
          n_proposed_swapE_moves <- n_proposed_swapE_moves + 1
          n_accepted_swapE_moves <- n_accepted_swapE_moves + tmp$accept
          if(tmp$accept==1)
          {
            curr_aug_dat <- tmp$new_aug_dat # if accepted move, update accordingly

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
    if(MCMC_settings$moves_switch$zeta_on)
    {
      tmp <- move_zeta_gibbs(curr_aug_dat,
                             curr_theta,
                             obs_dat,
                             hyperparameters)
      curr_theta <- tmp$new_theta # always update with new theta (Gibbs sampler)
    }

    # move mu
    #print("Move mu")
    if(MCMC_settings$moves_switch$mu_on)
    {
      for(g in seq_len(n_groups))
      {
        #print(g)
        for(j in seq(2,ncol(curr_aug_dat$D[[g]]),1))
        {
          #print(j)
          tmp <- move_lognormal(what="mu", g, j-1, MCMC_settings$moves_options$sdlog_mu[[g]][[j-1]],
                                curr_aug_dat,
                                curr_theta,
                                obs_dat,
                                hyperparameters,
                                index_dates)
          n_proposed_mu_moves[[g]][j-1] <- n_proposed_mu_moves[[g]][j-1] + 1
          n_accepted_mu_moves[[g]][j-1] <- n_accepted_mu_moves[[g]][j-1] + tmp$accept
          if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
        }
      }
    }

    # move CV
    #print("Move cv")
    if(MCMC_settings$moves_switch$CV_on)
    {
      for(g in seq_len(n_groups))
      {
        #print(g)
        for(j in seq(2,ncol(curr_aug_dat$D[[g]]),1))
        {
          #print(j)
          tmp <- move_lognormal(what="CV", g, j-1, MCMC_settings$moves_options$sdlog_CV[[g]][[j-1]],
                                curr_aug_dat,
                                curr_theta,
                                obs_dat,
                                hyperparameters,
                                index_dates)
          n_proposed_CV_moves[[g]][j-1] <- n_proposed_CV_moves[[g]][j-1] + 1
          n_accepted_CV_moves[[g]][j-1] <- n_accepted_CV_moves[[g]][j-1] + tmp$accept
          if (tmp$accept == 1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
        }
      }
    }

    # recording value of parameters and corresponding posterior after all moves
    if (output_stuff)
    {
      idx <- (k - MCMC_settings$chain_properties$burnin) / MCMC_settings$chain_properties$record_every+1
      theta_chain[[idx]] <- curr_theta
      aug_dat_chain[[idx]] <- curr_aug_dat
      logpost_chain[idx] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, hyperparameters, index_dates, range_dates) #### CONSIDER DOING THIS USING SAPPLY AFTER THE WHOLE THING
    }
  }

  ###############################################
  ### Compute acceptance probabilities ###
  ###############################################

  accept_prob <- list(
    D_moves=n_accepted_D_moves / n_proposed_D_moves,
    E_moves=n_accepted_E_moves / n_proposed_E_moves,
    mu_moves=lapply(seq_len(n_groups), function(g) n_accepted_mu_moves[[g]] / n_proposed_mu_moves[[g]]),
    CV_moves=lapply(seq_len(n_groups), function(g) n_accepted_CV_moves[[g]] / n_proposed_CV_moves[[g]]),
    zeta_moves=1)

###############################################
  ## Check for convergence
###############################################
  convergence <- check_convergence(theta_chain)
  ###############################################
  ### Return list of outputs of interest ###
  ###############################################

  res <- list(
    theta_chain = theta_chain, aug_dat_chain = aug_dat_chain,
    logpost_chain = logpost_chain, accept_prob = accept_prob,
    convergence = convergence
  )

  res
}

##' Splits a single MCMC chain into two
## to pretend we have run multiple chains
##' @details This is a utility function
##' and should not have to be called directly
##' by the user.
##' @param chain numeric vector, samples from MCMC chain.
##' @importFrom coda as.mcmc
##' @importFrom coda mcmc.list
##' @return \code{mcmc.list} from {coda}
##'
##' @export
split_chain_in_two <- function(chain) {
  len <- length(chain)
  if (len %% 2 != 0) {
    c1 <- as.mcmc(chain[1:((len + 1)/2)])
    c2 <- as.mcmc(chain[((len + 1)/2):len])
  } else {
    c1 <- as.mcmc(chain[1:(len/2)])
    c2 <- as.mcmc(chain[(len/2):len])
  }
  mcmc.list(c1, c2)
}
##' Check convergence of MCMC chains for
##' all parameters
##'
##'
##' @title Check convergence of MCMC chain
##' @param x MCMC chain, this is \code{theta_chain} from \code{RunMCMC}
##' @param threshold numeric, set to 1.1. This can be configured
##' but it is probably best not to change it. Gelman-Rubin diagnostic
##' should be approximately 1 at convergence.
##' @return a list with the same structure as a single
##' element of the input. Each element of th return value is a
##' logical indicating whether the chain has converged
##' @author Sangeeta Bhatia
##' @importFrom coda gelman.diag
##' @export
check_convergence <- function(x, threshold = 1.1) {
  ## check that the upper ci is
  ## smaller than the threshold. If below
  ## threshold, convergence is TRUE if not
  ## FALSE
  f <- function(gr) {
    unname(gr$psrf[ , "Upper C.I."] < threshold)
  }
  ## length of theta_chain is (n_iter - burnin)/record_every
  ## each element is a list with components mu and CV
  ## mu and CV are lists of length = number of groups
  ## Each element of mu (or CV) is a sample for the
  ## corresponding parameter for the group.
  ## To check convergence, we pull out the
  ## samples for each parameter from x and
  ## then compute the gelman.rubin diagnostic
  ## Return a list with the same structure
  ## as one element of x, but with
  ## logicals instead of samples from
  ## the parameter distribution
  out <- x[[1]]
  ## Check convergence of zeta
  zeta_chain <- vapply(x, function(y) y$zeta, numeric(1))
  zeta_chain <- split_chain_in_two(zeta_chain)
  gr_diag <- gelman.diag(zeta_chain, confidence = 0.95)
  out$zeta <- f(gr_diag)
  ngroups <- length(x[[1]]$mu)
  ## mu and CV, do one at a time
  for (group in seq_len(ngroups)) {
    ## Each group might have a different
    ## number of delays
    ndelays <- length(x[[group]]$mu)
    for (delay in seq_len(ndelays)) {
      chain <- vapply(
        x, function(y) y$mu[[group]][[delay]],
        numeric(1)
      )
      gr_diag <- gelman.diag(split_chain_in_two(chain), confidence = 0.95)
      out$mu[[group]][[delay]] <- f(gr_diag)

      ## Same for CV for this group and this delay
      chain <- vapply(
        x, function(y) y$CV[[group]][[delay]],
        numeric(1)
      )
      gr_diag <- gelman.diag(split_chain_in_two(chain), confidence = 0.95)
      out$CV[[group]][[delay]] <- f(gr_diag)
    }
  }
  out
}
#' Plots the MCMC chains of parameters
#'
#' @param MCMCres The output of function \code{\link{RunMCMC}}.
#' @param theta_true A list of parameters to which the output chains should be compared. If not \code{NULL}, this should contain:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups=length(MCMCres$aug_dat_chain[[1]]$D)}. Each element of \code{mu} should be a scalar or vector giving the mean delay(s) in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar or vector giving the coefficient o variation of the delay(s) in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' @return Nothing. Only performs a plot.
#' @import graphics
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
plot_parameter_chains <- function(MCMCres, theta_true=NULL)
{
  par(mfrow=c(2, 5),mar=c(5, 6, 1, 1))

  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol )
  n_groups <- length(n_dates)

  iterations <- seq_len(length(MCMCres$theta_chain))
  if(!is.null(theta_true)) x_coord_simul <- max(iterations)*1.05

  # looking at the logposterior chain
  plot(MCMCres$logpost_chain, type="l", xlab="Iterations", ylab="Log posterior")

  # looking at mean delay
  group_idx <- 1 ##########################
  j <- 1
  mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] )
  plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(non hospitalised-alive group)", ylim=c(0, 20))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$mu[[group_idx]][j])
  par(xpd=FALSE)

  legend("topright", "Onset-Report", lty=1)
  group_idx <- 2 ##########################
  j <- 1
  mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] )
  plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(non hospitalised-dead group)", ylim=c(0, 20))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$mu[[group_idx]][j], col=j)
  par(xpd=FALSE)
  for(j in seq(2,(n_dates[group_idx]-1),1))
  {
    mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] )
    lines(mu, col=j)
    par(xpd=TRUE)
    if(!is.null(theta_true)) points(x_coord_simul,theta_true$mu[[group_idx]][j], col=j)
    par(xpd=FALSE)
  }
  legend("topright", c("Onset-Death", "Onset-Report"), lty=1, col=seq_len(n_dates[group_idx]))
  group_idx <- 3 ##########################
  j <- 1
  mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] )
  plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(hospitalised-alive group)", ylim=c(0, 20))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$mu[[group_idx]][j], col=j)
  par(xpd=FALSE)
  for(j in seq(2,(n_dates[group_idx]-1), 1))
  {
    mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] )
    lines(mu, col=j)
    par(xpd=TRUE)
    if(!is.null(theta_true)) points(x_coord_simul,theta_true$mu[[group_idx]][j], col=j)
    par(xpd=FALSE)
  }
  legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), lty=1, col=seq_len(n_dates[group_idx]))
  group_idx <- 4 ##########################
  j <- 1
  mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] )
  plot(mu, type="l", xlab="Iterations", ylab="mean delays\n(hospitalised-dead group)", ylim=c(0, 20))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$mu[[group_idx]][j], col=j)
  par(xpd=FALSE)
  for(j in seq(2,(n_dates[group_idx]-1), 1))
  {
    mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] )
    lines(mu, col=j)
    par(xpd=TRUE)
    if(!is.null(theta_true)) points(x_coord_simul,theta_true$mu[[group_idx]][j], col=j)
    par(xpd=FALSE)
  }
  legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"), lty=1, col=seq_len(n_dates[group_idx]))

  # looking at zeta
  zeta <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta )
  plot(zeta, type="l", xlab="Iterations", ylab="zeta")
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$zeta)
  par(xpd=FALSE)

  # looking at CV delay
  group_idx <- 1 ##########################
  j <- 1
  CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] )
  plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(non hospitalised-alive group)", ylim=c(0, 2))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$CV[[group_idx]][j], col=j)
  par(xpd=FALSE)
  legend("topright", "Onset-Report", lty=1)
  group_idx <- 2 ##########################
  j <- 1
  CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] )
  plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(non hospitalised-dead group)", ylim=c(0, 2))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$CV[[group_idx]][j], col=j)
  par(xpd=FALSE)
  for(j in seq(2,(n_dates[group_idx]-1), 1))
  {
    CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] )
    lines(CV, col=j)
    par(xpd=TRUE)
    if(!is.null(theta_true)) points(x_coord_simul,theta_true$CV[[group_idx]][j], col=j)
    par(xpd=FALSE)
  }
  legend("topright", c("Onset-Death", "Onset-Report"), lty=1, col=seq_len(n_dates[group_idx]))
  group_idx <- 3 ##########################
  j <- 1
  CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] )
  plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(hospitalised-alive group)", ylim=c(0, 2))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$CV[[group_idx]][j], col=j)
  par(xpd=FALSE)
  for(j in seq(2,(n_dates[group_idx]-1),1))
  {
    CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] )
    lines(CV, col=j)
    par(xpd=TRUE)
    if(!is.null(theta_true)) points(x_coord_simul,theta_true$CV[[group_idx]][j], col=j)
    par(xpd=FALSE)
  }
  legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), lty=1, col=seq_len(n_dates[group_idx]))
  group_idx <- 4 ##########################
  j <- 1
  CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] )
  plot(CV, type="l", xlab="Iterations", ylab="CV delays\n(hospitalised-dead group)", ylim=c(0, 2))
  par(xpd=TRUE)
  if(!is.null(theta_true)) points(x_coord_simul,theta_true$CV[[group_idx]][j], col=j)
  par(xpd=FALSE)
  for(j in seq(2,(n_dates[group_idx]-1), 1))
  {
    CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] )
    lines(CV, col=j)
    par(xpd=TRUE)
    if(!is.null(theta_true)) points(x_coord_simul,theta_true$CV[[group_idx]][j], col=j)
    par(xpd=FALSE)
  }
  legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"), lty=1, col=seq_len(n_dates[group_idx]))
}


#' Plots the MCMC chains of augmented data
#'
#' @param MCMCres The output of function \code{\link{RunMCMC}}.
#' @param aug_dat_true A list containing the data to which the output chains should be compared. If not \code{NULL}, this should have the format of the output of function \code{\link{simul_true_data}}.
#' @return Nothing. Only performs a plot.
#' @import graphics
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
plot_aug_dat_chains <- function(MCMCres, aug_dat_true=NULL)
{
  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol )
  n_groups <- length(n_dates)
  n_indiv_per_group <- sapply(MCMCres$aug_dat_chain[[1]]$D, nrow )

  iterations <- seq_len(length(MCMCres$theta_chain))
  if(!is.null(aug_dat_true))
  {
    x_coord_simul <- max(iterations)*c(1.05, 1.07, 1.09)
    pch_types <- c(18, 21, 13)
    cex <- 1.5
  }

  par(mfrow=c(4, 5),mar=c(5, 6, 1, 1))
  group_idx <- 1 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for(i in seq_len(length(indiv_to_plot)) )
  {
    j <- 1
    date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
    plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
    par(xpd=TRUE)
    if(!is.null(aug_dat_true))
    {
      pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
      points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
    }
    par(xpd=FALSE)
    for(j in seq_len(n_dates[group_idx]))
    {
      date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
      lines(date, col=j)
      par(xpd=TRUE)
      if(!is.null(aug_dat_true))
      {
        pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
        points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
      }
      par(xpd=FALSE)
    }
    legend("topright", c("Onset","Report"), lty=1, col=seq_len(n_dates[group_idx]))
  }
  group_idx <- 2 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for(i in seq_len(length(indiv_to_plot)) )
  {
    j <- 1
    date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
    plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
    par(xpd=TRUE)
    if(!is.null(aug_dat_true))
    {
      pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
      points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
    }
    par(xpd=FALSE)
    for(j in seq(2, (n_dates[group_idx]), 1))
    {
      date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
      lines(date, col=j)
      par(xpd=TRUE)
      if(!is.null(aug_dat_true))
      {
        pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
        points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
      }
      par(xpd=FALSE)
    }
    legend("topright", c("Onset","Death","Report"), lty=1, col=seq_len(n_dates[group_idx]) )
  }
  group_idx <- 3 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for(i in seq_len(length(indiv_to_plot)) )
  {
    j <- 1
    date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
    plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
    par(xpd=TRUE)
    if(!is.null(aug_dat_true))
    {
      pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
      points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
    }
    par(xpd=FALSE)
    for(j in seq(2, (n_dates[group_idx]), 1))
    {
      date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
      lines(date, col=j)
      par(xpd=TRUE)
      if(!is.null(aug_dat_true))
      {
        pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
        points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
      }
      par(xpd=FALSE)
    }
    legend("topright", c("Onset","Hosp","Disch","Report"), lty=1, col=seq_len(n_dates[group_idx]) )
  }
  group_idx <- 4 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for(i in seq_len(length(indiv_to_plot)) )
  {
    j <- 1
    date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
    plot(date, type="l", xlab="Iterations", ylab="", ylim=c(min(date)-30, max(date)+30))
    par(xpd=TRUE)
    if(!is.null(aug_dat_true))
    {
      pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
      points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
    }
    par(xpd=FALSE)
    for(j in seq(2,(n_dates[group_idx]), 1))
    {
      date <- sapply(iterations, function(e) MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j] )
      lines(date, col=j)
      par(xpd=TRUE)
      if(!is.null(aug_dat_true))
      {
        pch <- pch_types[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))]
        points(x_coord_simul[match(aug_dat_true$E[[group_idx]][indiv_to_plot[i],j], c(-1, 0, 1))],aug_dat_true$D[[group_idx]][indiv_to_plot[i],j], col=j, pch=pch, cex=cex)
      }
      par(xpd=FALSE)
    }
    legend("topright", c("Onset","Hosp","Death","Report"), lty=1, col=seq_len(n_dates[group_idx]) )
  }
}

#' Compute correlation between the MCMC chains of mean and CV of each delay
#'
#' @param MCMCres The output of function \code{\link{RunMCMC}}.
#' @param plot A boolean indicating whether to plot the correlations or not
#' @return A list of results of correlation test (obtained from the function \code{\link{cor.test}}) between the posterior mean and the posterior CV for each delay.
#' @import graphics
#' @import stats
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
compute_correlations_mu_CV <- function(MCMCres, plot=TRUE)
{
  cor_mu_CV <- list()

  if(plot) par(mfrow=c(2, 5),mar=c(5, 6, 1, 1))

  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol )
  n_groups <- length(n_dates)

  iterations <- seq_len(length(MCMCres$theta_chain) )

  group_idx <- 1
  mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]])
  CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]])
  if(plot) plot(mu, CV, type="l")
  cor_mu_CV[[group_idx]] <- cor.test(mu, CV)

  for(group_idx in seq(2, n_groups, 1))
  {
    cor_mu_CV[[group_idx]] <- list()
    for(j in seq_len(n_dates[[group_idx]]-1))
    {
      mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
      CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
      if(plot) plot(mu, CV, type="l", col=j)
      cor_mu_CV[[group_idx]][[j]] <- cor.test(mu, CV)
    }
  }

  return(cor_mu_CV)
}

#' Compute autocorrelation for each parameter of the MCMC chains
#'
#' @param MCMCres The output of function \code{\link{RunMCMC}}.
#' @return A list of results of autocorrelation results (obtained from the function \code{\link{acf}}).
#' @import graphics
#' @import stats
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
compute_autocorr <- function(MCMCres)
{
  autocorr <- list()
  autocorr$mu <- list()
  autocorr$CV <- list()

  par(mfrow=c(4, 5),mar=c(4, 4, 4, 0.5))

  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol )
  n_groups <- length(n_dates)

  iterations <- seq_len(length(MCMCres$theta_chain) )

  for(group_idx in seq(1, n_groups, 1))
  {
    autocorr$mu[[group_idx]] <- list()
    autocorr$CV[[group_idx]] <- list()
    for(j in seq_len(n_dates[[group_idx]]-1))
    {
      mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
      CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
      autocorr$mu[[group_idx]][[j]] <- acf(mu, main=sprintf("Mu, group %d, delay %d", group_idx, j))
      autocorr$CV[[group_idx]][[j]] <- acf(CV, main=sprintf("CV, group %d, delay %d", group_idx, j))
    }
  }

  zeta <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta)
  autocorr$zeta <- acf(zeta, main="zeta")

  return(autocorr)
}


#' Computes posterior estimates of parameters from the MCMC chain
#'
#' @param MCMCres The output of function \code{\link{RunMCMC}}.
#' @param central A character specifying what the central estimate should be (median or mean posterior)
#' @param CrI A scalar in [0;1] used to compute the posterior credible intervals. For 95\% credible intervals, use CrI=0.95.
#' @param theta_true A list of parameters to which the output chains should be compared. If not \code{NULL}, this should contain:
#' \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups=length(MCMCres$aug_dat_chain[[1]]$D)}. Each element of \code{mu} should be a scalar or vector giving the mean delay(s) in that group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a scalar or vector giving the coefficient o variation of the delay(s) in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a data point is not missing, it is recorded with error.}
#' }
#' The posterior distributions of parameters are then plotted together with \code{theta_true}.
#' @param plot A boolean specifying whether to plot boxplots of the posterior estimates or not
#' @param cex.axis A numerical value giving the amount by which x axis labels should be magnified relative to the default.
#' @return A list containing two elements: the posterior estimates of parameters:
#' \itemize{
#'  \item{\code{logpost}}{: A vector of three values giving the central log-posterior estimate (first value) and quantiles corresponding to CrI (second and third values). }
#'  \item{\code{theta}}{: A list giving posterior parameter estimates
#'  \itemize{
#'  \item{\code{mu}}{: A list of length \code{n_groups=length(MCMCres$aug_dat_chain[[1]]$D)}. Each element of \code{mu} should be a matrix with 3 rows giving the posterior mean delay(s) in that group (1st row = central posterior estimate, 2nd and 3rd rows = credible interval) .}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of \code{CV} should be a matrix with 3 rows giving the posterior CV of the delay(s) in that group (1st row = central posterior estimate, 2nd and 3rd rows = credible interval) .}
#'  \item{\code{zeta}}{: A vector of three values in [0;1] giving the posterior estimate of the probability that, if a data point is not missing, it is recorded with error (1st value = central posterior estimate, 2nd and 3rd values = credible interval).}
#'  }
#'  }
#' }
#' @import graphics
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
get_param_posterior_estimates <- function(MCMCres, central=c("median","mean"), CrI=0.95, theta_true=NULL, plot=TRUE, cex.axis=1)
{
  par(mfrow=c(2, 5),mar=c(3, 5, 0.5, 0.5))

  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol )
  n_groups <- length(n_dates)
  n_indiv_per_group <- sapply(MCMCres$aug_dat_chain[[1]]$D, nrow )

  iterations <- seq_len(length(MCMCres$theta_chain))
  output <- list()

  # looking at the logposterior chain
  output$logpost <- c(get(central)(MCMCres$logpost_chain), quantile(MCMCres$logpost_chain, c((1-CrI)/2, CrI+(1-CrI)/2)) )
  if(plot)
  {
    boxplot(MCMCres$logpost_chain, ylab="Log Posterior", border="black", axes=FALSE)
    axis(side=1, at=1, labels="Log Posterior", tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
  }

  output$theta <- list()

  # looking at mean delay
  group_idx <- 1 ##########################
  mu <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] ) )
  output$theta$mu[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(mu[[j]]), quantile(mu[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(mu, ylab="mean delays\n(non hospitalised-alive group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels="Onset-Report", tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$mu[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }
  group_idx <- 2 ##########################
  mu <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] ) )
  output$theta$mu[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(mu[[j]]), quantile(mu[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(mu, ylab="mean delays\n(non hospitalised-dead group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels=c("Onset-Death", "Onset-Report"), tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$mu[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }
  group_idx <- 3 ##########################
  mu <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] ) )
  output$theta$mu[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(mu[[j]]), quantile(mu[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(mu, ylab="mean delays\n(hospitalised-alive group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels=c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$mu[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }
  group_idx <- 4 ##########################
  mu <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j] ) )
  output$theta$mu[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(mu[[j]]), quantile(mu[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(mu, ylab="mean delays\n(hospitalised-dead group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels=c("Onset-Hosp", "Hosp-Death", "Onset-Report"), tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$mu[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }

  # looking at zeta
  zeta <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta)
  if(plot)
  {
    boxplot(zeta, axes=FALSE, ylab="zeta")
    axis(side=1, at=1, labels="zeta", tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(theta_true$zeta, pch=8, lwd=2, cex=2)
  }

  # looking at CV delay
  group_idx <- 1 ##########################
  CV <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] ) )
  output$theta$CV[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(CV[[j]]), quantile(CV[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(CV, ylab="CV delays\n(non hospitalised-alive group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels="Onset-Report", tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$CV[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }
  group_idx <- 2 ##########################
  CV <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] ) )
  output$theta$CV[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(CV[[j]]), quantile(CV[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(CV, ylab="CV delays\n(non hospitalised-dead group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels=c("Onset-Death", "Onset-Report"), tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$CV[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }
  group_idx <- 3 ##########################
  CV <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] ) )
  output$theta$CV[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(CV[[j]]), quantile(CV[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(CV, ylab="CV delays\n(hospitalised-alive group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels=c("Onset-Hosp", "Hosp-Disch", "Onset-Report"), tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$CV[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }
  group_idx <- 4 ##########################
  CV <- lapply(seq_len(n_dates[group_idx]-1), function(j) sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j] ) )
  output$theta$CV[[group_idx]] <- sapply(seq_len(n_dates[group_idx]-1), function(j) c(get(central)(CV[[j]]), quantile(CV[[j]], c((1-CrI)/2, CrI+(1-CrI)/2)) ) )
  if(plot)
  {
    boxplot(CV, ylab="CV delays\n(hospitalised-dead group)", main="", border=seq_len(n_dates[group_idx]-1), axes=FALSE)
    axis(side=1, at=seq_len(n_dates[group_idx]-1), labels=c("Onset-Hosp", "Hosp-Death", "Onset-Report"), tick=FALSE, cex.axis=cex.axis)
    axis(side=2)
    if(!is.null(theta_true)) points(seq_len(n_dates[group_idx]-1), theta_true$CV[[group_idx]], pch=8, lwd=2, cex=2, col=seq_len(n_dates[group_idx]-1))
  }

  return(output)

}

