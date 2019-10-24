###############################################
### functions to handle dates ###
###############################################

#' Convert date to integer corresponding to number of days from a given origin
#'
#' @param date \code{Date} object to be converted to integer.
#' @param origin The date used as origin for counting days.
#' @return An integer corresponding to the number of days since \code{origin}.
#' @export
#' @examples
#' date_to_int(as.Date("1970-01-01"))
date_to_int <- function(date, origin = "1970-01-01")
{
  return(as.integer(date - as.Date(origin)))
}

#' Convert integer to date, based on a given origin from which the integer counts the number of days
#' 
#' @param int integer to be converted to \code{Date}.
#' @param origin The date used as origin for counting days.
#' @return A \code{Date} object corresponding to \code{int} days after \code{origin}.
#' @export
#' @examples
#' int_to_date(365)
int_to_date <- function(int, origin = "1970-01-01")
{
  return(int + as.Date(origin))
}

###############################################
### functions to find parameters of distributions from mean/var or mean/sd ###
###############################################

#' Finds the two shape parameters of the beta distribution, given a mean and variance
#' 
#' @param mean Mean of the beta distribution
#' @param var Variance of the beta distribution
#' @return A vector containing the two shape parameters of the beta distribution.
#' @export
#' @examples
#' param_beta <- find_params_beta(0.1, 0.05)
#' sample_beta <- rbeta(1000, param_beta[1], param_beta[2])
#' mean(sample_beta) # compare to 0.1
#' var(sample_beta) # compare to 0.05
find_params_beta <- function(mean, var) # function to determine parameters of the beta distribution corresponding to a given mean and variance
{
  # for a beta distribution:
  # mean = shape1/(shape1+shape2) 
  # var = shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1)) 
  # after solving we find that
  shape1 <- mean^2*(1-mean)/var - mean
  shape2 <- shape1*(1/mean-1)
  return(c(shape1, shape2))
}

#' Finds the shape and scale parameters of the gamma distribution, given a mean and standard deviation
#' 
#' @param mean Mean of the gamma distribution
#' @param sigma Standard deviation of the gamma distribution
#' @param CV An alternative way to specify the stndard deviation, through the coefficient of variation, that is \code{sigma}/\code{mean}
#' @return A vector containing the shape and scale parameters of the gamma distribution.
#' @export
#' @examples
#' param_gamma <- find_params_gamma(mean=0.1, sigma=0.05)
#' sample_gamma <- rgamma(1000, shape=param_gamma[1], scale=param_gamma[2])
#' mean(sample_gamma) # compare to 0.1
#' sd(sample_gamma) # compare to 0.05
find_params_gamma <- function(mean, sigma=mean*CV, CV) # function to determine parameters of the gamma distribution corresponding to a given mean and std
{
  shape <- (mean/sigma)^2
  scale <- sigma^2/(mean)
  return(c(shape, scale))
}

###############################################
### compute_delta functions to compute relevant delays based on index, which tells you which dates should be used for delayl calculation ###
###############################################

#' Compute relevant delays for one individual, based on data
#' 
#' @param D A list of data, in the format of the first element (called \code{true_dat}) in the list returned by \code{\link{simul_true_data}}. 
#' @param group_idx a scalar or vector giving the index of the group(s) for which to calculate the delays
#' @param indiv_idx a scalar or vector giving the index of the individuals in group \code{group_idx} for which to calculate the delays
#' @param delay_idx a scalar or vector giving the index of the delay(s) to be calculated
#' @param index_dates A list containing indications on which delays to consider in the simulation. 
#' @details \code{index_dates} should be a list of length \code{n_groups=length(D)}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. 
#' For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' @return A list of same length as \code{D}. Each element in the list is a matrix with same number of rows as in \code{D}, but potentially a different number of columns, corresponding to the relevant delays calculated for that group. 
#' @export
#' @examples
#' ### Number of groups of individuals to simulate ###
#' n_groups <- 2
#' ### Number of dates to simulate for each group ###
#' n_dates <- c(2, 3)
#' ### Setting up the parameters for the simulation ###
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of these delays
#' ### Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each group? ###
#' index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)))
#' ### Perform the simulation ###
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' ### Compute the first delay for first individual in first group ###
#' compute_delta_group_delay_and_indiv(D$true_dat, group_idx=1, indiv_idx=1, delay_idx=1, index_dates)
compute_delta_group_delay_and_indiv<-function(D, group_idx, indiv_idx, delay_idx, index_dates)
{
  Delta <- D[[group_idx]][indiv_idx,index_dates[[group_idx]][,delay_idx][2]] - D[[group_idx]][indiv_idx,index_dates[[group_idx]][,delay_idx][1]]
  return(Delta)
}

#' Compute relevant delays based on data
#' 
#' @param D A list of data, in the format of the first element (called \code{true_dat}) in the list returned by \code{\link{simul_true_data}}. 
#' @param index_dates A list containing indications on which delays to consider in the simulation. 
#' @details \code{index_dates} should be a list of length \code{n_groups=length(D)}. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. 
#' For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' @return A list of same length as \code{D}. Each element in the list is a matrix with same number of rows as in \code{D}, but potentially a different number of columns, corresponding to the relevant delays calculated for that group. 
#' @export
#' @examples
#' ### Number of groups of individuals to simulate ###
#' n_groups <- 2
#' ### Number of dates to simulate for each group ###
#' n_dates <- c(2, 3)
#' ### Setting up the parameters for the simulation ###
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of these delays
#' ### Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each group? ###
#' index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)))
#' ### Perform the simulation ###
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' ### Compute the corresponding delays ###
#' Delays <- compute_delta(D$true_dat, index_dates)
compute_delta <- function(D, index_dates)
{
  Delta <- lapply(seq_len(length(D)), function(g){
    m <- matrix(NA, nrow(D[[g]]), ncol(D[[g]])-1)
    for(j in seq_len( ncol(m)) )
    {
      m[,j] <- D[[g]][,index_dates[[g]][,j][2]] - D[[g]][,index_dates[[g]][,j][1]]
    }
    return(m)
  })
  return(Delta)
}

###############################################
### find range of dates
###############################################

#' Find range of dates from a dataset
#' 
#' @param obs_dat A list of data, in the format of the first element (called \code{obs_dat}) in the list returned by \code{\link{simul_obs_dat}}. 
#' @return A vector of two integers coresponding to the range in \code{obs_dat}. 
#' @import stats
#' @export
#' @examples
#' ### Number of groups of individuals to simulate ###
#' n_groups <- 2
#' ### Number of dates to simulate for each group ###
#' n_dates <- c(2, 3)
#' ### Setting up the parameters for the simulation ###
#' theta <- list()
#' theta$mu <- list(5, c(10, 15)) # mean delays, for each group
#' theta$CV <- list(0.5, c(0.5, 0.5)) # coefficient of variation of these delays
#' ### Number of individuals to simulate in each group ###
#' n_per_group <- rep(10, n_groups)
#' ### Range of dates in which to draw the first set of dates for each group ###
#' range_dates <- date_to_int(c(as.Date("01/01/2014", "%d/%m/%Y"), as.Date("01/01/2015", "%d/%m/%Y")))
#' ### Which delays to use to simulate subsequent dates from the first, in each group? ###
#' index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)))
#' ### Perform the simulation ###
#' D <- simul_true_data(theta, n_per_group, range_dates, index_dates)
#' ### Find the range ###
#' find_range(D$true_dat)
#' ### Compare with range specified in the first place ### 
#' range_dates
find_range <- function(obs_dat)
{
  min_date <- min(obs_dat[[1]][,1], na.rm=TRUE)
  max_date <- max(obs_dat[[1]][,1], na.rm=TRUE)
  for(g in seq_len(length(obs_dat)) )
  {
    for(j in seq_len(ncol(obs_dat[[g]])) )
    {
      min_date_tmp <- min(obs_dat[[g]][,j], na.rm=TRUE)
      min_date <- min(c(min_date, min_date_tmp), na.rm=TRUE)
      
      max_date_tmp <- max(obs_dat[[g]][,j], na.rm=TRUE)
      max_date <- max(c(max_date, max_date_tmp), na.rm=TRUE)
    }
  }
  return(c(min_date, max_date))
}


###############################################
### compute rules on the order of dates based on index_dates object
###############################################

#' Compute rules on the order of dates based on index_dates object; see details
#' 
#' @param index_dates A list containing indications on which delays to consider in the simulation, see details.
#' @details \code{index_dates} should be a list; each elements corresponding to a group of individuals of interest. Each element of \code{index_dates} should be a matrix with 2 rows and a number of columns corresponding to the delays of interest for that group. For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date. 
#' 
#' If index_dates[[k]] has two columns containing respectively c(1, 2) and c(1, 3), this indicates that for group \code{k} we are interested in two delays: the first delay being between date 1 and date 2, and the second being between date 1 and date 3. 
#' 
#' This function is used to find appropriate starting points for the MCMC, i.e. when choosing initial values for the missing dates, this allows making sure the chosen value is consistent with the ordering of dates in that group, and hence will not generate a null likelihood. 
#' 
#' @return A list of same lenght as index_dates, containing indications on ordering of dates for each group. More specifically, each element of \code{index_dates_order} is a matrix with 2 rows and a number of columns corresponding to the delays with order rules for that group. 
#' For each column (i.e. each delay), the first row gives the index of the origin date, and the second row gives the index of the destination date.
#' Each column specifies a rule saying that the origin date must be before the destination date.  
#' This is computed from \code{index_dates} assuming all delays specified in \code{index_dates} have to be positive, and using transitivity rules to derive potential additional rules of positivity. 
#' In the example below, in group 3, \code{index_dates} indicates that the delay between dates 1 and 2, and the delay between dates 2 and 3, are positive. Hence by transitivity, the delay between date 1 and 3 has to be positive as well, as seen in the output of the function in that example. 
#' @export
#' @examples
#' index_dates <- list(matrix(c(1, 2), nrow=2), 
#'                        cbind(c(1, 2), c(1, 3)), 
#'                        cbind(c(1, 2), c(2, 3), c(1, 4)), 
#'                        cbind(c(1, 2), c(2, 3), c(1, 4)) )
#' index_dates_order <- compute_index_dates_order(index_dates)
compute_index_dates_order <- function(index_dates)
{
  index_dates_order <- index_dates
  
  number_cols_added_at_this_round <- 0
  for(e in 1:length(index_dates))
  {
    tmp <- index_dates_order[[e]]
    link <- tmp[2,tmp[2,] %in% tmp[1,]]
    for(k in link)
    {
      for(i in which(tmp[2,]==k))
      {
        for(j in which(tmp[1,]==k))
        {
          if(!any(sapply(1:ncol(tmp), function(e) all(tmp[,e] == c(tmp[1,i],tmp[2,j])) ))) # this means the rule obtained by transitivity is not yet present --> needs to be added
          {
            index_dates_order[[e]] <- cbind(index_dates_order[[e]], c(tmp[1,i],tmp[2,j]) )
            number_cols_added_at_this_round <- number_cols_added_at_this_round + 1
          }
        }
      }
    }
  }
  while(number_cols_added_at_this_round>0)
  {
    number_cols_added_at_this_round <- 0
    for(e in 1:length(index_dates))
    {
      tmp <- index_dates_order[[e]]
      link <- tmp[2,tmp[2,] %in% tmp[1,]]
      for(k in link)
      {
        for(i in which(tmp[2,]==k))
        {
          for(j in which(tmp[1,]==k))
          {
            if(!any(sapply(1:ncol(tmp), function(e) all(tmp[,e] == c(tmp[1,i],tmp[2,j])) ))) # this means the rule obtained by transitivity is not yet present --> needs to be added
            {
              index_dates_order[[e]] <- cbind(index_dates_order[[e]], c(tmp[1,i],tmp[2,j]) )
              number_cols_added_at_this_round <- number_cols_added_at_this_round + 1
            }
          }
        }
      }
    }
  }
  return(index_dates_order)
}






check_MCMC_settings <- function(MCMC_settings, index_dates)
{
  if(MCMC_settings$chain_properties$n_iter < MCMC_settings$chain_properties$burnin)
    stop("Burnin must be <= n_iter")
  if(length(MCMC_settings$moves_options$sdlog_mu) != length(index_dates))
    stop("sdlog_mu does not have the correct length")
  if(length(MCMC_settings$moves_options$sdlog_CV) != length(index_dates))
    stop("sdlog_CV does not have the correct length")
  if(!all(lengths(MCMC_settings$moves_options$sdlog_mu) == lengths(index_dates) / 2))
    stop("sdlog_mu does not have the correct structure")
  if(!all(lengths(MCMC_settings$moves_options$sdlog_CV) == lengths(index_dates) / 2))
    stop("sdlog_CV does not have the correct structure")
  if(!is.numeric(MCMC_settings$tol) | (MCMC_settings$tol<0) | (MCMC_settings$tol>1))
    stop("tol should be a number between 0 and 1")
}


infer_missing_dates <- function(D, 
                                E = NULL, 
                                g, # index of the group
                                e, # index of individual in that group
                                index_dates_order, 
                                do_not_infer_from = NULL, 
                                theta = NULL, 
                                tol = 1e-6, 
                                move_to_D = NULL) # used to compute probabilities of a given move; this only contains the values of the dates for that specific individual
{
  #if(g == 3 & e == 71)
  #{
  #  browser()
  #}
  
  D_proxy <- D
  if(is.null(E))
  {
    true_missing_dates <- which(is.na(D_proxy[[g]][e,]))
  } else 
  {
    #browser()
    true_missing_dates <- which(E[[g]][e,] == -1)
  }
  missing_dates <- union(true_missing_dates, do_not_infer_from)
  other_missing_dates <- setdiff(missing_dates, true_missing_dates)
  D_proxy[[g]][e,missing_dates] <- NA
  
  #print(true_missing_dates)
  #print(other_missing_dates)
  
  if(is.null(theta))
  {
    while(length(true_missing_dates)>0)
    {
      ### basic inference using only order of dates
      can_be_inferred_from <- lapply(missing_dates, function(i) {
        x <- which(index_dates_order[[g]]==i, arr.ind = TRUE)
        from_idx <- sapply(seq_len(nrow(x)), function(k) index_dates_order[[g]][-x[k,1],x[k,2]] )
        from_value <- sapply(seq_len(nrow(x)), function(k) D_proxy[[g]][e,index_dates_order[[g]][-x[k,1],x[k,2]]])
        rule <- sapply(seq_len(nrow(x)), function(k) if(x[k,1]==1) "before" else "after"  )
        return(list(rule=rule,from_idx=from_idx, from_value=from_value))
      })
      idx_can_be_inferred <- which(sapply(seq_len(length(missing_dates)), 
                                          function(i) any(!is.na(can_be_inferred_from[[i]]$from_value))))
      for(k in idx_can_be_inferred)
      {
        x <- which(!is.na(can_be_inferred_from[[k]]$from_value))
        if(length(x)==1)
        {
          inferred <- can_be_inferred_from[[k]]$from_value[x]
        }else
        {
          if(all(can_be_inferred_from[[k]]$rule[x] == "before"))
          {
            inferred <- min(can_be_inferred_from[[k]]$from_value[x])
          }else if (all(can_be_inferred_from[[k]]$rule[x] == "after"))
          {
            inferred <- max(can_be_inferred_from[[k]]$from_value[x])
          }else
          {
            max_val <- min(can_be_inferred_from[[k]]$from_value[x][can_be_inferred_from[[k]]$rule[x] %in% "before"])
            min_val <- max(can_be_inferred_from[[k]]$from_value[x][can_be_inferred_from[[k]]$rule[x] %in% "after"])
            if(min_val>max_val)
            {
              stop("Incompatible data to infer from. ")
            }else
            {
              inferred <- floor(median(c(min_val, max_val)))
            }
          }
        }
        D_proxy[[g]][e,missing_dates[k]] <- inferred
      }
      missing_dates <- which(is.na(D_proxy[[g]][e,]))
      true_missing_dates <- setdiff(missing_dates, other_missing_dates)
      #print("rep")
      #print(missing_dates)
      #print(true_missing_dates)
      #print(other_missing_dates)
      #Sys.sleep(0.1)
    }
    prob <- 1
    if(!is.null(move_to_D))
    {
      to_compare <- which(!is.na(D_proxy[[g]][e,]))
      if(all(D_proxy[[g]][e,to_compare] == move_to_D[to_compare])) 
      {
        prob_move_to_D <- 1 
      } else 
      {
        prob_move_to_D <- 0
      }
    } else
    {
      prob_move_to_D <- NA
    }
  }else
  {
    prob <- rep(NA, length(D_proxy[[g]][e,]))
    prob_move_to_D <- rep(NA, length(D_proxy[[g]][e,]))
    while(length(true_missing_dates)>0)
    {
      ### more sophisticated inference using the delays
      can_be_inferred_directly_from <- lapply(missing_dates, function(date_idx) {
        infer_directly_from(g, date_idx, D_proxy, e, theta) })
      idx_can_be_inferred_directly <- which(sapply(seq_len(length(missing_dates)), 
                                                   function(i) length(can_be_inferred_directly_from[[i]]$from_value)>0))
      ### TO DO: even more sophisticated allowing for sums of delays if needed
      for(k in idx_can_be_inferred_directly)
      {
        # browser()
        tmp <- sample_new_date_value(g, can_be_inferred_directly_from[[k]], index_dates, tol = tol)
        D_proxy[[g]][e,missing_dates[k]] <- tmp$inferred
        prob[missing_dates[k]] <- tmp$probability_inferred_value
        if(!is.null(move_to_D)) 
        {
          if(any(tmp$all_possible_values == move_to_D[missing_dates[k]]))
          {
            prob_move_to_D[missing_dates[k]] <- tmp$probabilities[tmp$all_possible_values == move_to_D[missing_dates[k]]]
          } else
          {
            prob_move_to_D[missing_dates[k]] <- 0
          }
        }
      }
      missing_dates <- which(is.na(D_proxy[[g]][e,]))
      true_missing_dates <- setdiff(missing_dates, other_missing_dates)
      #print("rep")
      #print(missing_dates)
      #print(true_missing_dates)
      #print(other_missing_dates)
      #Sys.sleep(0.1)
    }
    
  }
  
  to_replace <- is.na(D_proxy[[g]][e, ])
  if(any(to_replace)) D_proxy[[g]][e, to_replace] <- D[[g]][e, to_replace]
  prob <- prod(prob[!is.na(prob)])
  if(all(is.na(prob_move_to_D)))
  {
    prob_move_to_D <- NA 
  } else
  {
    prob_move_to_D <- prod(prob_move_to_D[!is.na(prob_move_to_D)])
  }
  
  res <- list(D = D_proxy, 
              prob = prob, 
              prob_move_to_D = prob_move_to_D)
  
  return(res)
}

# for a given group, date and individual, this proposes a new value for that date. 
# NOTE this assumes that there are parameters available to infer delays from
# needs to first get can_be_inferred_directly_from from infer_directly_from function
# The function returns a list containing:
# - inferred: the proposed value
# - probability_inferred_value the probability this specific value was proposed (rather than other values)
# all_possible_values: all values which could have been proposed,
# probabilities: probabilities of all values which could have been proposed
sample_new_date_value <- function(g, # group index
                                  can_be_inferred_directly_from, # the output of function infer_directly_from
                                  index_dates, 
                                  tol = 1e-6) # used to define the tail of the CDF of the delay - anything in the tail is neglected
{
  can_be_inferred_directly_from <- can_be_inferred_directly_from[which(!is.na(can_be_inferred_directly_from$from_value)), ]
  
  tmp <- t(sapply(1:nrow(can_be_inferred_directly_from), 
                  function(ii) find_params_gamma(
                    can_be_inferred_directly_from$mu[ii], 
                    CV = can_be_inferred_directly_from$CV[ii])))
  can_be_inferred_directly_from$shape_delay <- tmp[,1]
  can_be_inferred_directly_from$scale_delay <- tmp[,2]
  delay_max <- max(
    ceiling(qgamma(1 - tol, 
                   shape = can_be_inferred_directly_from$shape_delay, 
                   scale = can_be_inferred_directly_from$scale_delay)))
  
  prob_delay <- sapply(1:nrow(can_be_inferred_directly_from), 
                       function(ii) pgamma(0:(delay_max+1), 
                                           shape = can_be_inferred_directly_from$shape_delay[ii], 
                                           scale = can_be_inferred_directly_from$scale_delay[ii]))
  weights <- apply(prob_delay, 2, diff)
  delays <- 0:(delay_max)
  
  all_possible_values <- matrix(rep(can_be_inferred_directly_from$from_value, length(delays)), nrow = length(delays), byrow = TRUE) + 
    matrix(rep(can_be_inferred_directly_from$multiply, length(delays)), nrow = length(delays), byrow = TRUE) * 
    matrix(rep(delays, nrow(can_be_inferred_directly_from)), nrow = length(delays), byrow = FALSE)
  
  if(nrow(can_be_inferred_directly_from) > 1)
  {
    #need to merge info from various delays
    
    merged_possible_values <- all_possible_values[,1]
    for(jj in 2:ncol(all_possible_values))
    {
      merged_possible_values <- intersect(merged_possible_values, all_possible_values[,jj])
    }
    
    if(length(merged_possible_values)==0)
    {
      warning("Incompatible dates, drawing from the first delay only.")
      return(sample_new_date_value(g, can_be_inferred_directly_from[1,], 
                                   index_dates, tol) 
      )
    } 
    
    merged_weights <- sapply(merged_possible_values, function(ee) prod(weights[all_possible_values == ee]))
    if(sum(merged_weights)==0)
    {
      warning("Incompatible dates, drawing from the first delay only.")
      return(sample_new_date_value(g, can_be_inferred_directly_from[1,], 
                                   index_dates, tol) 
      )
    } 
    
    merged_weights <- merged_weights / sum(merged_weights)
    
    tmp <- rmultinom(1, 1, merged_weights)
    
    inferred <- sum(merged_possible_values*tmp)
    
    probability_inferred_value <- sum(merged_weights*tmp) ### needs to be used for calculating the probability of proposing each value?
    
    return(list(inferred = inferred,
                probability_inferred_value = probability_inferred_value,
                all_possible_values = merged_possible_values,
                probabilities = merged_weights))
    
  }else
  {
    #draw randomly from the only possible delay 
    tmp <- sapply(1:nrow(can_be_inferred_directly_from), 
                  function(ii) rmultinom(1, 1, weights[,ii]))
    
    inferred <- colSums(all_possible_values*tmp)
    
    probability_inferred_value <- colSums(weights * tmp) ### needs to be used for calculating the probability of proposing each value?
    
    return(list(inferred = inferred,
                probability_inferred_value = probability_inferred_value,
                all_possible_values = all_possible_values,
                probabilities = weights))
    
  }
}

# for a given group, date and individual, returns a list showing where that 
# specific date can be inferred from. The list contains:
# - multiply (+ or -1 depending on whether the date to infer is after or before 
# the date(s) it can be inferred from)
# - from_idx: index of the date(s) it can be inferred from
# - from_value: current value of the date(s) it can be inferred from
# - mu: mean delay comparef to the date(s) it can be inferred from
# - CV: CV delay comparef to the date(s) it can be inferred from
infer_directly_from <- function(g, # group
                                date_idx, # date index for that group
                                D, # data (augmented)
                                i, # index of individual in the group
                                theta) # parameter values 
{
  delay_idx <- which(sapply(1:ncol(index_dates[[g]]), function(k) date_idx %in% index_dates[[g]][,k]))
  from_idx <- sapply(delay_idx, function(k) index_dates[[g]][,k][index_dates[[g]][,k] != date_idx])
  from_value <- D[[g]][i, from_idx]
  if(!is.null(theta) )
  {
    mu <- theta$mu[[g]][delay_idx] 
    CV <- theta$CV[[g]][delay_idx]
  } else 
  {
    mu <- NULL
    CV <- NULL
  }
  multiply <- sapply(delay_idx, function(e)
  {
    if(date_idx == index_dates[[g]][,e][1]) return(-1) else return(1)
  })
  to_keep <- which(!is.na(from_value))
  from_idx <- from_idx[to_keep]
  from_value <- from_value[to_keep]
  mu <- mu[to_keep]
  CV <- CV[to_keep]
  multiply <- multiply[to_keep]
  return(data.frame(multiply=multiply,
                    from_idx=from_idx, 
                    from_value=from_value,
                    mu = mu,
                    CV = CV))
}

my_mode <- function(x) { # function to get the mode
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

get_consensus <- function(aug_dat_true, MCMCres, posterior = c("mode", "median"))
{
  posterior <- match.arg(posterior)
  consensus_D <- aug_dat_true$D
  consensus_E <- aug_dat_true$E
  for(g in 1:length(consensus_D))
  {
    for(i in 1:nrow(consensus_D[[g]]))
    {
      for(j in 1:ncol(consensus_D[[g]]))
      {
        if(posterior == "median")
        {
          consensus_D[[g]][i,j] <- median(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, j]))
          consensus_E[[g]][i,j] <- median(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][i, j]))
        } else if(posterior == "mode")
        {
          consensus_D[[g]][i,j] <- my_mode(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$D[[g]][i, j]))
          consensus_E[[g]][i,j] <- my_mode(sapply(seq_len(length(MCMCres$aug_dat_chain)), function(e) MCMCres$aug_dat_chain[[e]]$E[[g]][i, j]))
        }
      }
    }
  }
  return(list(D = consensus_D, E = consensus_E))
}

### Compute sensitivity and specificity of detecting errors in dates
# remove missing dates from the denominator
compute_sensitivity_specificity_from_consensus <- function(aug_dat_true, consensus)
{
  tmp <- aug_dat_true$E
  aug_dat_true_E_no_missing <- tmp
  consensus$E_no_missing <- consensus$E
  for(g in 1:length(tmp))
  {
    aug_dat_true_E_no_missing[[g]] [tmp[[g]] == -1] <- NA 
    consensus$E_no_missing[[g]] [consensus$E[[g]] == -1] <- NA
  }
  
  sensitivity <- specificity <- rep(NA, 4)
  false_pos <- false_neg <- list()
  
  for(g in 1:length(tmp))
  {
    tab <- table(aug_dat_true_E_no_missing[[g]], consensus$E_no_missing[[g]])
    n_true_neg <- tab["0", "0"]
    n_true_pos <- tab["1", "1"]
    n_false_pos <- tab["0", "1"]
    n_false_neg <- tab["1", "0"]
    
    sensitivity[g] <- n_true_pos / (n_true_pos + n_false_neg)
    specificity[g] <- n_true_neg / (n_true_neg + n_false_pos)
    
    false_pos[[g]] <- which((aug_dat_true_E_no_missing[[g]] == 0) & # are really NOT an error
                              (consensus$E_no_missing[[g]]  == 1), # and are detected as errors
                            arr.ind = TRUE)
    
    false_neg[[g]] <- which((aug_dat_true_E_no_missing[[g]] == 1) & # are really  an error
                              (consensus$E_no_missing[[g]]  == 0), # and are NOT detected as errors
                            arr.ind = TRUE)
  }
  
  return(list(sensitivity = sensitivity,
              specificity = specificity,
              false_pos = false_pos,
              false_neg = false_neg))
}
