###############################################
### define augmented data to be used for initialisation of the chain ###
###############################################

### D contains the unobserved true dates ###

#' @import stats
initialise_aug_data <- function(obs_dat, index_dates_order)
{
  n_groups <- length(obs_dat)
  D <- list()
  for(g in 1:n_groups) 
  {
    D[[g]] <- obs_dat[[g]]
    for(e in 1:nrow(D[[g]]))
    {
      #print(e)
      
      # first deal with incompatible dates
      for(j in 1:ncol(index_dates_order[[g]]))
      {
        if(!any(is.na(D[[g]][e,index_dates_order[[g]][,j]])))
        {
          if(D[[g]][e,index_dates_order[[g]][1,j]] > D[[g]][e,index_dates_order[[g]][2,j]]) # there is a problem if the dates are in the wrong order
          {
            # check if there is one of the dates involved in more than one problematic delays, if so must be the problematic one:
            tmp <- table(as.vector(index_dates_order[[g]][,sapply(1:ncol(index_dates_order[[g]]), function(j) D[[g]][e,index_dates_order[[g]][1,j]] > D[[g]][e,index_dates_order[[g]][2,j]] )]))
            if(any(tmp>1))
            {
              must_be_wrong <- which.max(tmp)[1]
            }else
            {
              # check which of all dates is most outlier compared to all other dates, and if several take the first one as the wrong one
              diff_from_median <- abs(D[[g]][e,index_dates_order[[g]][,j]] - median(D[[g]][e,], na.rm=TRUE))
              must_be_wrong <- which(diff_from_median %in% max(diff_from_median))[1]
              must_be_wrong <- index_dates_order[[g]][,j][must_be_wrong]
            }
            D[[g]][e,must_be_wrong] <- NA
            while(!(must_be_wrong %in% index_dates_order[[g]][,j]))
            {
              # check if there is one of the dates involved in more than one problematic delays, if so must be the problematic one:
              tmp <- table(as.vector(index_dates_order[[g]][,sapply(1:ncol(index_dates_order[[g]]), function(j) D[[g]][e,index_dates_order[[g]][1,j]] > D[[g]][e,index_dates_order[[g]][2,j]] )]))
              if(any(tmp>1))
              {
                must_be_wrong <- which.max(tmp)[1]
              }else
              {
                # check which of all dates is most outlier compared to all other dates, and if several take the first one as the wrong one
                diff_from_median <- abs(D[[g]][e,index_dates_order[[g]][,j]] - median(D[[g]][e,], na.rm=TRUE))
                must_be_wrong <- which(diff_from_median %in% max(diff_from_median))[1]
                must_be_wrong <- index_dates_order[[g]][,j][must_be_wrong]
              }
              D[[g]][e,must_be_wrong] <- NA
            }
          }
        }
      }
      
      # now deal with missing dates
      missing_dates <- which(is.na(D[[g]][e,]))
      while(length(missing_dates)>0)
      {
        can_be_inferred_from <- lapply(missing_dates, function(i) {
          x <- which(index_dates_order[[g]]==i, arr.ind = TRUE)
          from_idx <- sapply(1:nrow(x), function(k) index_dates_order[[g]][-x[k,1],x[k,2]] )
          from_value <- sapply(1:nrow(x), function(k) D[[g]][e,index_dates_order[[g]][-x[k,1],x[k,2]]])
          rule <- sapply(1:nrow(x), function(k) if(x[k,1]==1) "before" else "after"  )
          return(list(rule=rule,from_idx=from_idx, from_value=from_value))
        })
        can_be_inferred <- which(sapply(1:length(missing_dates), function(i) any(!is.na(can_be_inferred_from[[i]]$from_value))))
        for(k in can_be_inferred)
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
          D[[g]][e,missing_dates[k]] <- inferred
        }
        missing_dates <- which(is.na(D[[g]][e,]))
      }
    }
  }
  names(D) <- names(obs_dat)
  
  # compute E accordingly
  E <- list()
  for(g in 1:n_groups) 
  {
    E[[g]] <- as.data.frame(matrix(NA,nrow(obs_dat[[g]]),ncol(obs_dat[[g]])))
    for(j in 1:ncol(obs_dat[[g]]))
    {
      E[[g]][!(D[[g]][,j] %in% obs_dat[[g]][,j]),j] <- 1 # error
      E[[g]][D[[g]][,j] %in% obs_dat[[g]][,j],j] <- 0 # no error
      E[[g]][is.na(obs_dat[[g]][,j]),j] <- -1 # missing value
    }
    names(E[[g]]) <- names(obs_dat[[g]])
  }
  names(E) <- names(obs_dat)
  
  aug_dat <- list(D = D,
                  E = E)
  
  return(aug_dat)
}

###############################################
### define parameters to be used for initialisation of the chain ###
###############################################

initialise_theta_from_aug_dat <- function(aug_dat, index_dates, zeta_init=0.1) # zeta_init doesn't really matter as we then use Gibbs sampler so will move fast to better values
{
  n_groups <- length(aug_dat$D)
  n_dates <- sapply(aug_dat$D, ncol)
    
  ### mean and std of distribution of various delays, by group
  ### we use a the starting point the observed mean and std of each delay in each group
  obs_delta <- compute_delta(aug_dat$D, index_dates)
  mu <- lapply(1:n_groups, function(g) abs(apply(obs_delta[[g]], 2, mean, na.rm=TRUE) ))
  names(mu) <- names(n_dates)
  sigma <- lapply(1:n_groups, function(g) abs(apply(obs_delta[[g]], 2, sd, na.rm=TRUE) ))
  names(sigma) <- names(n_dates)
  CV <- lapply(1:n_groups, function(g) sigma[[g]]/mu[[g]])
  names(CV) <- names(n_dates)
  
  ### list of all parameters
  theta <- list(zeta = zeta_init, # zeta is the probability for a date to be misrecorded, conditional on being recorded (<-> Ei != - 1)
                # TODO:
                # could consider having zeta being type of date specific (e.g. more error on onset than death dates),
                # time specific and/or space specific
                mu = mu, # mean of gamma distributions used to characterise the various delays in different groups: mu[[g]][k] is the mean k^th delay in group g
                CV = CV) # CV of gamma distributions used to characterise the various delays in different groups: CV[[g]][k] is the CV k^th delay in group g
  
  return(theta)
  
}
