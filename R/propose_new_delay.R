choose_delay <- function(i,
                              group_idx,
                              date_idx,
                              curr_aug_dat,
                              theta,
                              obs_dat,
                              hyperparameters,
                              index_dates,
                              range_dates = NULL) {

  # Identify delays this date is involved in ----------------------------------

  # Identify where the date is in index_dates (returns row and col numbers)
  x <- which(index_dates[[group_idx]] == date_idx, arr.ind = TRUE)

  # Identify delay(s) that this particular date is involved in
  which_delay <- x[, 2]

  # For each relevant delay, find the other date it's paired with
  from_idx <- sapply(seq_len(nrow(x)),
                     function(k) index_dates[[group_idx]][-x[k, 1], x[k, 2]])

  # Extract the values of those paired dates for individual i
  from_value <- sapply(seq_len(nrow(x)),
                       function(k) curr_aug_dat$D[[group_idx]][i, from_idx[k]])

  # If several delays involve this date, choose one at random
  if (length(from_idx) > 1) {
    tmp <- sample(seq_len(length(from_idx)), 1)
    which_delay <- which_delay[tmp]
    from_idx <- from_idx[tmp]
    from_value <- if (length(i) > 1) from_value[, tmp] else from_value[tmp]
  }

  if (date_idx < from_idx) {
    curr_delay <- from_value - curr_aug_dat$D[[group_idx]][i, date_idx]
  } else {
    curr_delay <- curr_aug_dat$D[[group_idx]][i, date_idx] - from_value
  }

  return(list(which_delay = which_delay, from_idx = from_idx, from_value = from_value, curr_delay = curr_delay))

}

propose_new_delay <- function(i,
                              group_idx,
                              date_idx,
                              curr_aug_dat,
                              theta,
                              obs_dat,
                              hyperparameters,
                              index_dates,
                              range_dates = NULL) {

  # Identify delays this date is involved in ----------------------------------

  tmp <- choose_delay(i, group_idx, date_idx, curr_aug_dat, theta, obs_dat,
                           hyperparameters, index_dates, range_dates)
  which_delay <- tmp$which_delay
  from_idx <- tmp$from_idx
  from_value <- tmp$from_value
  curr_delay <- tmp$curr_delay

  # Sample a new delay using the discretised Gamma distribution ---------------
  sample_delay <- discr_gamma_sample(length(i),
                                     mu = theta$mu[[group_idx]][which_delay],
                                     cv = theta$CV[[group_idx]][which_delay])

  # Depending on whether this date is before of after its pair, add or subtract
  ## ANNE: I think we do this multiple times across the code so would be good to wrap in single function
  if (date_idx < from_idx) {
    proposed_aug_dat_value <- from_value - sample_delay
  } else {
    proposed_aug_dat_value <- from_value + sample_delay
  }

  ## ANNE: need to select one option:
  # option 1 is keep doing this meaning we can propose the same value as before and this can be accepted / rejected and it's all the same
  # option 2 use a while loop to make sure the new value is different, but then we may need to do sth more clever for the acceptance probability

  # Create a copy of augmented data and insert the proposed value
  proposed_aug_dat <- curr_aug_dat
  proposed_aug_dat$D[[group_idx]][i, date_idx] <- proposed_aug_dat_value

  # Update error indicators accordingly ---------------------------------------
  # i.e. if D_i moves to y_i then E_i moves to 0, else E_i moves to 1.

  # Identify missing (if y_i missing then E_i = -1)
  missing <- which(is.na(obs_dat[[group_idx]][i, date_idx]))
  proposed_aug_dat$E[[group_idx]][i, date_idx][missing] <- -1

  # Identify non-erroneous (y_i observed without error then E_i = 0)
  non_erroneous <- which(proposed_aug_dat$D[[group_idx]][i, date_idx] ==
                           obs_dat[[group_idx]][i, date_idx])
  proposed_aug_dat$E[[group_idx]][i, date_idx][non_erroneous] <- 0

  # Identify erroneous observations (y_i observed with error then E_i = 1)
  erroneous <- which(!is.na(obs_dat[[group_idx]][i, date_idx]) &
                       proposed_aug_dat$D[[group_idx]][i, date_idx] !=
                       obs_dat[[group_idx]][i, date_idx])
  proposed_aug_dat$E[[group_idx]][i, date_idx][erroneous] <- 1

  return(list(proposed_aug_dat = proposed_aug_dat,
              curr_delay = curr_delay,
              sample_delay = sample_delay,
              which_delay = which_delay))
}



get_correct_factor_new_delay <- function(curr_delay, new_delay, theta, group_idx, which_delay)
{

  prob_proposing_curr_value <- DiscrGamma(curr_delay,
                                          mu = theta$mu[[group_idx]][which_delay],
                                          cv = theta$CV[[group_idx]][which_delay],
                                          log = TRUE)

  prob_proposing_new_value <- DiscrGamma(new_delay,
                                         mu = theta$mu[[group_idx]][which_delay],
                                         cv = theta$CV[[group_idx]][which_delay],
                                         log = TRUE)

  return(c(prob_proposing_curr_value = prob_proposing_curr_value,
           prob_proposing_new_value = prob_proposing_new_value))
}
