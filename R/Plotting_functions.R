#' Plots the MCMC chains of parameters 
#' 
#' @param MCMCres The output of function \code{\link{RunMCMC}}. 
#' @param theta_true A list of parameters to which the output chains should be
#'  compared. If not \code{NULL}, this should contain:
#' \itemize{
#'  \item{\code{mu}}{: A list of length
#'   \code{n_groups = length(MCMCres$aug_dat_chain[[1]]$D)}. Each element of
#'    \code{mu} should be a scalar or vector giving the mean delay(s) in that
#'     group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar or vector giving the coefficient of variation
#'    of the delay(s) in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' @return Nothing. Only performs a plot.
#' @import graphics
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
plot_parameter_chains <- function(MCMCres, theta_true = NULL) {
  
  par(mfrow = c(2, 5), mar = c(5, 6, 1, 1))
  
  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol)
  n_groups <- length(n_dates)
  
  iterations <- seq_len(length(MCMCres$theta_chain))
  if (!is.null(theta_true)) x_coord_simul <- max(iterations) * 1.05
  
  # looking at the logposterior chain 
  plot(MCMCres$logpost_chain, type = "l",
       xlab = "Iterations", ylab = "Log posterior")
  
  # looking at mean delay
  
  group_idx <- 1 ##########################
  j <- 1
  mu <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
  
  plot(mu, type = "l",
       xlab = "Iterations",
       ylab = "mean delays\n(non hospitalised-alive group)",
       ylim = c(0, 20))
  
  par(xpd = TRUE)
  
  if (!is.null(theta_true)) points(x_coord_simul, theta_true$mu[[group_idx]][j])
  
  par(xpd = FALSE)
  
  legend("topright", "Onset-Report", lty = 1)
  
  group_idx <- 2 ##########################
  j <- 1
  mu <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
  plot(mu, type = "l",
       xlab = "Iterations",
       ylab = "mean delays\n(non hospitalised-dead group)",
       ylim = c(0, 20))
  
  par(xpd = TRUE)
  
  if (!is.null(theta_true)) {
    points(x_coord_simul, theta_true$mu[[group_idx]][j], col = j)
  }
  
  par(xpd = FALSE)
  
  for (j in seq(2, (n_dates[group_idx] - 1), 1)) {
    
    mu <- sapply(iterations,
                 function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
    lines(mu, col = j)
    par(xpd = TRUE)
    
    if (!is.null(theta_true)) {
      points(x_coord_simul, theta_true$mu[[group_idx]][j], col = j)
    }
    
    par(xpd = FALSE)
  }
  
  legend("topright", c("Onset-Death", "Onset-Report"),
         lty = 1, col = seq_len(n_dates[group_idx]))
  
  group_idx <- 3 ##########################
  j <- 1
  mu <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
  
  plot(mu, type = "l",
       xlab = "Iterations",
       ylab = "mean delays\n(hospitalised-alive group)",
       ylim = c(0, 20))
  
  par(xpd = TRUE)
  if (!is.null(theta_true)) {
    points(x_coord_simul, theta_true$mu[[group_idx]][j], col = j)
  }
  
  par(xpd = FALSE)
  
  for (j in seq(2, (n_dates[group_idx] - 1), 1)) {
    mu <- sapply(iterations,
                 function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
    lines(mu, col = j)
    par(xpd = TRUE)
    if (!is.null(theta_true)) {
      points(x_coord_simul, theta_true$mu[[group_idx]][j], col = j)
    }
    par(xpd = FALSE)
  }
  legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"),
         lty = 1, col = seq_len(n_dates[group_idx]))
  
  group_idx <- 4 ##########################
  j <- 1
  mu <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
  
  plot(mu, type = "l", 
       xlab = "Iterations",
       ylab = "mean delays\n(hospitalised-dead group)",
       ylim = c(0, 20))
  
  par(xpd = TRUE)
  
  if (!is.null(theta_true)) {
    points(x_coord_simul, theta_true$mu[[group_idx]][j], col = j)
  }
  
  par(xpd = FALSE)
  
  for (j in seq(2, (n_dates[group_idx] - 1), 1)) {
    
    mu <- sapply(iterations,
                 function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
    lines(mu, col = j)
    par(xpd = TRUE)
    if (!is.null(theta_true)) {
      points(x_coord_simul, theta_true$mu[[group_idx]][j], col = j) 
    }
    par(xpd = FALSE)
  }
  
  legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"),
         lty = 1, col = seq_len(n_dates[group_idx]))
  
  # looking at zeta
  zeta <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta)
  plot(zeta, type = "l", xlab = "Iterations", ylab = "zeta")
  par(xpd = TRUE)
  if (!is.null(theta_true)) points(x_coord_simul, theta_true$zeta)
  par(xpd = FALSE)
  
  # looking at CV delay
  group_idx <- 1 ##########################
  j <- 1
  CV <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
  plot(CV, type = "l",
       xlab = "Iterations",
       ylab = "CV delays\n(non hospitalised-alive group)",
       ylim = c(0, 2))
  par(xpd = TRUE)
  if (!is.null(theta_true)) {
    points(x_coord_simul, theta_true$CV[[group_idx]][j], col = j)
  }
  par(xpd = FALSE)
  legend("topright", "Onset-Report", lty = 1)
  
  group_idx <- 2 ##########################
  j <- 1
  CV <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
  plot(CV, type = "l",
       xlab = "Iterations",
       ylab = "CV delays\n(non hospitalised-dead group)",
       ylim = c(0, 2))
  par(xpd = TRUE)
  if (!is.null(theta_true)) {
    points(x_coord_simul, theta_true$CV[[group_idx]][j], col = j)
  }
  par(xpd = FALSE)
  for (j in seq(2,(n_dates[group_idx] - 1), 1)) {
    CV <- sapply(iterations,
                 function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
    lines(CV, col = j)
    par(xpd = TRUE)
    if (!is.null(theta_true)) {
      points(x_coord_simul, theta_true$CV[[group_idx]][j], col = j)
    }
    par(xpd = FALSE)
  }
  legend("topright", c("Onset-Death", "Onset-Report"),
         lty = 1, col = seq_len(n_dates[group_idx]))
  
  group_idx <- 3 ##########################
  j <- 1
  CV <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
  plot(CV, type = "l",
       xlab = "Iterations",
       ylab = "CV delays\n(hospitalised-alive group)",
       ylim = c(0, 2))
  par(xpd = TRUE)
  if (!is.null(theta_true)) {
    points(x_coord_simul, theta_true$CV[[group_idx]][j], col = j)
  }
  par(xpd = FALSE)
  for (j in seq(2, (n_dates[group_idx] - 1), 1)) {
    CV <- sapply(iterations,
                 function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
    lines(CV, col = j)
    par(xpd = TRUE)
    if (!is.null(theta_true)) {
      points(x_coord_simul, theta_true$CV[[group_idx]][j], col = j)
    }
    par(xpd = FALSE)
  }
  legend("topright", c("Onset-Hosp", "Hosp-Disch", "Onset-Report"),
         lty = 1, col = seq_len(n_dates[group_idx]))
  
  group_idx <- 4 ##########################
  j <- 1
  CV <- sapply(iterations,
               function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
  plot(CV, type = "l",
       xlab = "Iterations",
       ylab = "CV delays\n(hospitalised-dead group)",
       ylim = c(0, 2))
  par(xpd = TRUE)
  if(!is.null(theta_true)) {
    points(x_coord_simul, theta_true$CV[[group_idx]][j], col = j)
  }
  par(xpd = FALSE)
  for(j in seq(2, (n_dates[group_idx] - 1), 1)) {
    CV <- sapply(iterations,
                 function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
    lines(CV, col = j)
    par(xpd = TRUE)
    if (!is.null(theta_true)) {
      points(x_coord_simul, theta_true$CV[[group_idx]][j], col = j)
    }
    par(xpd = FALSE)
  }
  legend("topright", c("Onset-Hosp", "Hosp-Death", "Onset-Report"),
         lty = 1, col = seq_len(n_dates[group_idx]))
}


#' Plots the MCMC chains of augmented data 
#' 
#' @param MCMCres The output of function \code{\link{RunMCMC}}. 
#' @param aug_dat_true A list containing the data to which the output chains
#'  should be compared. If not \code{NULL}, this should have the format of the
#'   output of function \code{\link{simul_true_data}}. 
#' @return Nothing. Only performs a plot.
#' @import graphics
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
plot_aug_dat_chains <- function(MCMCres, aug_dat_true = NULL) {
  
  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol)
  n_groups <- length(n_dates)
  n_indiv_per_group <- sapply(MCMCres$aug_dat_chain[[1]]$D, nrow)
  
  iterations <- seq_len(length(MCMCres$theta_chain))
  if (!is.null(aug_dat_true)) {
    x_coord_simul <- max(iterations) * c(1.05, 1.07, 1.09)
    pch_types <- c(18, 21, 13)
    cex <- 1.5
  }
  
  par(mfrow = c(4, 5), mar = c(5, 6, 1, 1))
  
  group_idx <- 1 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for (i in seq_len(length(indiv_to_plot))) {
    j <- 1
    date <- sapply(iterations, function(e) {
      MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
    })
    plot(date, type = "l",
         xlab = "Iterations",
         ylab = "",
         ylim = c(min(date) - 30, max(date) + 30))
    par(xpd = TRUE)
    if (!is.null(aug_dat_true)) {
      pch <- pch_types[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )]
      points(x_coord_simul[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )],
      aug_dat_true$D[[group_idx]][indiv_to_plot[i], j], 
      col = j, pch = pch, cex = cex)
    }
    par(xpd = FALSE)
    for (j in seq_len(n_dates[group_idx])) {
      date <- sapply(iterations, function(e) {
        MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
      })
      lines(date, col = j)
      par(xpd = TRUE)
      if (!is.null(aug_dat_true)) {
        pch <- pch_types[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )]
        points(x_coord_simul[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )],
        aug_dat_true$D[[group_idx]][indiv_to_plot[i], j],
        col = j, pch = pch, cex = cex)
      }
      par(xpd = FALSE)
    }
    legend("topright", c("Onset", "Report"),
           lty = 1, col = seq_len(n_dates[group_idx]))
  }
  
  group_idx <- 2 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for (i in seq_len(length(indiv_to_plot))) {
    j <- 1
    date <- sapply(iterations, function(e) {
      MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
    })
    plot(date, type = "l",
         xlab = "Iterations",
         ylab = "",
         ylim = c(min(date) - 30, max(date) + 30))
    par(xpd = TRUE)
    if (!is.null(aug_dat_true)) {
      pch <- pch_types[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )]
      points(x_coord_simul[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )],
      aug_dat_true$D[[group_idx]][indiv_to_plot[i], j],
      col = j, pch = pch, cex = cex)
    }
    par(xpd = FALSE)
    for (j in seq(2, (n_dates[group_idx]), 1)) {
      date <- sapply(iterations, function(e) {
        MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
      })
      lines(date, col = j)
      par(xpd = TRUE)
      if (!is.null(aug_dat_true)) {
        pch <- pch_types[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )]
        points(x_coord_simul[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )],
        aug_dat_true$D[[group_idx]][indiv_to_plot[i], j],
        col = j, pch = pch, cex = cex)
      }
      par(xpd = FALSE)
    }
    legend("topright", c("Onset", "Death", "Report"),
           lty = 1, col = seq_len(n_dates[group_idx]))
  }
  
  group_idx <- 3 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for (i in seq_len(length(indiv_to_plot))) {
    j <- 1
    date <- sapply(iterations, function(e) {
      MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
    })
    plot(date, type = "l",
         xlab = "Iterations",
         ylab = "",
         ylim = c(min(date) - 30, max(date) + 30))
    par(xpd = TRUE)
    if (!is.null(aug_dat_true)) {
      pch <- pch_types[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )]
      points(x_coord_simul[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )],
      aug_dat_true$D[[group_idx]][indiv_to_plot[i], j],
      col = j, pch = pch, cex = cex)
    }
    par(xpd = FALSE)
    for (j in seq(2, (n_dates[group_idx]), 1)) {
      date <- sapply(iterations, function(e) {
        MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
      })
      lines(date, col = j)
      par(xpd = TRUE)
      if (!is.null(aug_dat_true)) {
        pch <- pch_types[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )]
        points(x_coord_simul[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )],
        aug_dat_true$D[[group_idx]][indiv_to_plot[i], j],
        col = j, pch = pch, cex = cex)
      }
      par(xpd = FALSE)
    }
    legend("topright", c("Onset", "Hosp", "Disch", "Report"),
           lty = 1, col = seq_len(n_dates[group_idx]))
  }
  
  group_idx <- 4 ##########################
  # randomly pick 5 individuals in that group
  indiv_to_plot <- sample(seq_len(n_indiv_per_group[group_idx]), 5)
  for (i in seq_len(length(indiv_to_plot))) {
    j <- 1
    date <- sapply(iterations, function(e) {
      MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
    })
    plot(date, type = "l",
         xlab = "Iterations",
         ylab = "",
         ylim = c(min(date) - 30, max(date) + 30))
    par(xpd = TRUE)
    if (!is.null(aug_dat_true)) {
      pch <- pch_types[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )]
      points(x_coord_simul[match(
        aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
      )],
      aug_dat_true$D[[group_idx]][indiv_to_plot[i], j],
      col = j, pch = pch, cex = cex)
    }
    par(xpd = FALSE)
    for (j in seq(2, (n_dates[group_idx]), 1)) {
      date <- sapply(iterations, function(e) {
        MCMCres$aug_dat_chain[[e]]$D[[group_idx]][indiv_to_plot[i], j]
      })
      lines(date, col = j)
      par(xpd = TRUE)
      if (!is.null(aug_dat_true)) {
        pch <- pch_types[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )]
        points(x_coord_simul[match(
          aug_dat_true$E[[group_idx]][indiv_to_plot[i], j], c(-1, 0, 1)
        )],
        aug_dat_true$D[[group_idx]][indiv_to_plot[i], j],
        col = j, pch = pch, cex = cex)
      }
      par(xpd = FALSE)
    }
    legend("topright", c("Onset", "Hosp", "Death", "Report"),
           lty = 1, col = seq_len(n_dates[group_idx]))
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
compute_correlations_mu_CV <- function(MCMCres, plot = TRUE) {
  
  cor_mu_CV <- list()
  
  if (plot) par(mfrow = c(2, 5), mar = c(5, 6, 1, 1))
  
  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol)
  n_groups <- length(n_dates)
  
  iterations <- seq_len(length(MCMCres$theta_chain))
  
  group_idx <- 1
  mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]])
  CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]])
  if (plot) plot(mu, CV, type = "l")
  cor_mu_CV[[group_idx]] <- cor.test(mu, CV)
  
  for (group_idx in seq(2, n_groups, 1)) {
    cor_mu_CV[[group_idx]] <- list()
    for (j in seq_len(n_dates[[group_idx]] - 1)) {
      mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
      CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
      if (plot) plot(mu, CV, type = "l", col = j)
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
compute_autocorr <- function(MCMCres) {
  autocorr <- list()
  autocorr$mu <- list()
  autocorr$CV <- list()
  
  par(mfrow = c(4, 5), mar = c(4, 4, 4, 0.5))
  
  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol)
  n_groups <- length(n_dates)
  
  iterations <- seq_len(length(MCMCres$theta_chain))
  
  for (group_idx in seq(1, n_groups, 1)) {
    autocorr$mu[[group_idx]] <- list()
    autocorr$CV[[group_idx]] <- list()
    for (j in seq_len(n_dates[[group_idx]] - 1)) {
      mu <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
      CV <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
      autocorr$mu[[group_idx]][[j]] <- acf(mu,
                                           main = sprintf(
                                             "Mu, group %d, delay %d",
                                             group_idx,
                                             j))
      autocorr$CV[[group_idx]][[j]] <- acf(CV,
                                           main = sprintf(
                                             "CV, group %d, delay %d",
                                             group_idx,
                                             j))
    }
  }
  
  zeta <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta)
  autocorr$zeta <- acf(zeta, main = "zeta")
  
  return(autocorr)
}


#' Computes posterior estimates of parameters from the MCMC chain
#' 
#' @param MCMCres The output of function \code{\link{RunMCMC}}. 
#' @param central A character specifying what the central estimate should be
#'  (median or mean posterior)
#' @param CrI A scalar in [0;1] used to compute the posterior credible
#'  intervals. For 95\% credible intervals, use CrI=0.95.
#' @param theta_true A list of parameters to which the output chains should be
#'  compared. If not \code{NULL}, this should contain:
#' \itemize{
#'  \item{\code{mu}}{: A list of length
#'   \code{n_groups = length(MCMCres$aug_dat_chain[[1]]$D)}. Each element of
#'    \code{mu} should be a scalar or vector giving the mean delay(s) in that
#'     group.}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a scalar or vector giving the coefficient o variation
#'    of the delay(s) in that group.}
#'  \item{\code{zeta}}{: A scalar in [0;1] giving the probability that, if a
#'   data point is not missing, it is recorded with error.}
#' }
#' The posterior distributions of parameters are then plotted together with
#'  \code{theta_true}. 
#' @param plot A boolean specifying whether to plot boxplots of the posterior
#'  estimates or not
#' @param cex.axis A numerical value giving the amount by which x axis labels
#'  should be magnified relative to the default.
#' @return A list containing two elements: the posterior estimates of
#'  parameters:
#' \itemize{
#'  \item{\code{logpost}}{: A vector of three values giving the central
#'   log-posterior estimate (first value) and quantiles corresponding to CrI
#'    (second and third values).}
#'  \item{\code{theta}}{: A list giving posterior parameter estimates
#'  \itemize{
#'  \item{\code{mu}}{: A list of length
#'   \code{n_groups = length(MCMCres$aug_dat_chain[[1]]$D)}. Each element of
#'    \code{mu} should be a matrix with 3 rows giving the posterior mean
#'     delay(s) in that group (1st row = central posterior estimate, 2nd and
#'      3rd rows = credible interval).}
#'  \item{\code{CV}}{: A list of length \code{n_groups}. Each element of
#'   \code{CV} should be a matrix with 3 rows giving the posterior CV of the
#'    delay(s) in that group (1st row = central posterior estimate, 2nd and 3rd
#'     rows = credible interval) .}
#'  \item{\code{zeta}}{: A vector of three values in [0;1] giving the posterior
#'   estimate of the probability that, if a data point is not missing, it is
#'    recorded with error (1st value = central posterior estimate, 2nd and 3rd
#'     values = credible interval).}
#'  }
#'  }
#' }
#' @import graphics
#' @export
#' @examples
#' ### TO WRITE OR ALTERNATIVELY REFER TO VIGNETTE TO BE WRITTEN ###
get_param_posterior_estimates <- function(MCMCres,
                                          central = c("median", "mean"),
                                          CrI = 0.95,
                                          theta_true = NULL,
                                          plot = TRUE,
                                          cex.axis = 1) {
  
  par(mfrow = c(2, 5), mar = c(3, 5, 0.5, 0.5))
  
  n_dates <- sapply(MCMCres$aug_dat_chain[[1]]$D, ncol)
  n_groups <- length(n_dates)
  n_indiv_per_group <- sapply(MCMCres$aug_dat_chain[[1]]$D, nrow)
  
  iterations <- seq_len(length(MCMCres$theta_chain))
  output <- list()
  
  # looking at the logposterior chain 
  output$logpost <- c(
    get(central)(MCMCres$logpost_chain),
    quantile(MCMCres$logpost_chain, c((1 - CrI) / 2, CrI + (1 - CrI) / 2))
  )
  
  if (plot) {
    boxplot(MCMCres$logpost_chain,
            ylab = "Log Posterior",
            border = "black",
            axes = FALSE)
    axis(side = 1,
         at = 1,
         labels = "Log Posterior",
         tick = FALSE,
         cex.axis = cex.axis)
    axis(side = 2)
  }
  
  output$theta <- list()
  
  # looking at mean delay 
  group_idx <- 1 ##########################
  mu <- lapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
    }
  )
  output$theta$mu[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(mu[[j]]),
        quantile(mu[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    }
  )
  if (plot) {
    boxplot(mu,
            ylab = "mean delays\n(non hospitalised-alive group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1, at = seq_len(n_dates[group_idx] - 1),
         labels = "Onset-Report", tick = FALSE, cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$mu[[group_idx]],
             pch = 8, lwd = 2, cex = 2, col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  group_idx <- 2 ##########################
  mu <- lapply(
    seq_len(n_dates[group_idx] - 1), function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
    }
  )
  output$theta$mu[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(mu[[j]]),
        quantile(mu[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    }
  )
  if (plot) {
    boxplot(mu,
            ylab = "mean delays\n(non hospitalised-dead group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1, at = seq_len(n_dates[group_idx]-1),
         labels = c("Onset-Death", "Onset-Report"),
         tick = FALSE, cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$mu[[group_idx]],
             pch = 8, lwd = 2, cex = 2, col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  group_idx <- 3 ##########################
  mu <- lapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
    }
  )
  output$theta$mu[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(mu[[j]]),
        quantile(mu[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    }
  )
  if (plot) {
    boxplot(mu,
            ylab = "mean delays\n(hospitalised-alive group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1,
         at = seq_len(n_dates[group_idx] - 1),
         labels = c("Onset-Hosp", "Hosp-Disch", "Onset-Report"),
         tick = FALSE, cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$mu[[group_idx]],
             pch = 8, lwd = 2, cex = 2,
             col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  group_idx <- 4 ##########################
  mu <- lapply(
    seq_len(n_dates[group_idx]-1),
    function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$mu[[group_idx]][j])
    })
  output$theta$mu[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(mu[[j]]),
        quantile(mu[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    })
  if (plot) {
    boxplot(mu,
            ylab = "mean delays\n(hospitalised-dead group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1,
         at = seq_len(n_dates[group_idx] - 1),
         labels = c("Onset-Hosp", "Hosp-Death", "Onset-Report"),
         tick = FALSE,
         cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$mu[[group_idx]],
             pch = 8, lwd = 2, cex = 2,
             col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  # looking at zeta
  zeta <- sapply(iterations, function(e) MCMCres$theta_chain[[e]]$zeta)
  if (plot) {
    boxplot(zeta, axes = FALSE, ylab = "zeta")
    axis(side = 1, at = 1, labels = "zeta", tick = FALSE, cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) points(theta_true$zeta, pch = 8, lwd = 2, cex = 2)
  }
  
  # looking at CV delay 
  group_idx <- 1 ##########################
  CV <- lapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
    })
  output$theta$CV[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(CV[[j]]),
        quantile(CV[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    })
  if (plot) {
    boxplot(CV,
            ylab = "CV delays\n(non hospitalised-alive group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1,
         at = seq_len(n_dates[group_idx] - 1),
         labels = "Onset-Report",
         tick = FALSE,
         cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$CV[[group_idx]],
             pch = 8, lwd = 2, cex = 2, col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  group_idx <- 2 ##########################
  CV <- lapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
    })
  output$theta$CV[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(CV[[j]]),
        quantile(CV[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    })
  if (plot) {
    boxplot(CV,
            ylab = "CV delays\n(non hospitalised-dead group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1,
         at = seq_len(n_dates[group_idx] - 1),
         labels = c("Onset-Death", "Onset-Report"),
         tick = FALSE,
         cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$CV[[group_idx]],
             pch = 8, lwd = 2, cex = 2,
             col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  group_idx <- 3 ##########################
  CV <- lapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
    })
  output$theta$CV[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(CV[[j]]),
        quantile(CV[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    })
  if (plot) {
    boxplot(CV,
            ylab = "CV delays\n(hospitalised-alive group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1,
         at = seq_len(n_dates[group_idx] - 1),
         labels = c("Onset-Hosp", "Hosp-Disch", "Onset-Report"),
         tick = FALSE,
         cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$CV[[group_idx]],
             pch = 8, lwd = 2, cex = 2,
             col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  group_idx <- 4 ##########################
  CV <- lapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      sapply(iterations,
             function(e) MCMCres$theta_chain[[e]]$CV[[group_idx]][j])
    })
  output$theta$CV[[group_idx]] <- sapply(
    seq_len(n_dates[group_idx] - 1),
    function(j) {
      c(get(central)(CV[[j]]),
        quantile(CV[[j]], c((1 - CrI) / 2, CrI + (1 - CrI) / 2)))
    })
  if (plot) {
    boxplot(CV,
            ylab = "CV delays\n(hospitalised-dead group)",
            main = "",
            border = seq_len(n_dates[group_idx] - 1),
            axes = FALSE)
    axis(side = 1,
         at = seq_len(n_dates[group_idx] - 1),
         labels = c("Onset-Hosp", "Hosp-Death", "Onset-Report"),
         tick = FALSE,
         cex.axis = cex.axis)
    axis(side = 2)
    if (!is.null(theta_true)) {
      points(seq_len(n_dates[group_idx] - 1),
             theta_true$CV[[group_idx]],
             pch = 8, lwd = 2, cex = 2,
             col = seq_len(n_dates[group_idx] - 1))
    }
  }
  
  return(output)
  
}

