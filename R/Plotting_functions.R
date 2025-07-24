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

