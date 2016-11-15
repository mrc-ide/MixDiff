load("~/Dropbox/Ebola/DATA (WHO, population, demo)/Ebola Data/WHO Ebola Line List_2015-09-28/dbCombined_cleanv42_2015-09-28.Rdata")

### dates we're interested in
colDates <- c("DateOnset", "DateHospitalCurrentAdmit", "DateDeath", "DateDischargeHosp", "DateReport")

### keep only confirmed cases ###
dat <- dat[dat$EpiCaseDef %in% 1,] ##  is Confirmed

### remove those with missing hospitalisation status of final status ###
dat <- dat[(!is.na(dat$HospitalizedEver)) & !(is.na(dat$FinalStatus)),]

### checking that individuals who have a recorded date of hospitalisation or discharge are recorded as hospitalised
table(dat$HospitalizedEver,!is.na(dat$DateHospitalCurrentAdmit) | !is.na(dat$DateDischargeHosp) )

### checking that individuals who have a recorded date of death have final outcome death
table(dat$FinalStatus,!is.na(dat$DateDeath) )
dat$FinalStatus[!is.na(dat$DateDeath)] <- "Dead"

### classification of people according to their hospitalisation status and final outcome
table(dat$HospitalizedEver, dat$FinalStatus, useNA="ifany")
dat$Path <- paste0("Hosp",dat$HospitalizedEver,"-",dat$FinalStatus)

### dates we are interested in ###
# date of exposure ### could use dat$ContactDateStart1 if only one exposure reported and ContactDateStart1 = ContactDateEnd1??? OR COULD INFER DATE OF EXPOSURE AS WELL AMONGST SEVERAL / UNCERTAIN DATES
head(dat$DateOnset) # date of onset
head(dat$DateHospitalCurrentAdmit) # date of hospitalisation
head(dat$DateDeath) # date of death
head(dat$DateDischargeHosp) # date of discharge
head(dat$DateReport) # date of report

### creating the dataset to use ###
dat <- dat[,c("GlobalRecordId","DateOnset","DateHospitalCurrentAdmit","DateDeath","DateDischargeHosp","DateReport","Path")]

### remove individuals with no dates at all among the dates we're interested in ###
dat <- dat[!(rowSums(!is.na(dat[,colDates])) == 0),]

### making sure that individuals who are hospitalised and die have NA in their discharge date (not interested in the delay from discharge to death which doesn't really make sense) ###
dat$DateDischargeHosp[dat$Path == "HospYes-Dead"] <- NA # note we should never really looked at these anyway

### make sure date columns are in date format ###
for(e in colDates) dat[,e] <- as.Date(dat[,e])

###############################################
###############################################
### parameter estimation using MCMC ###
###############################################
###############################################

###############################################
### data ###
###############################################

dat<-readRDS("~/Desktop/Dat.rds")

# head(dat)

N <- nrow(dat) # Number of cases


tmp <- split(dat[colDates],dat$Path)
# splitting dataset according to Path and removing NA date columns in each of these
# - should only remain dates that are relevant for each group
dat_by_group <- lapply(tmp, function(x) x[,colSums(is.na(x))!=nrow(x)] )

n_dates <- sapply(dat_by_group, ncol )

###############################################
### define parameters ###
###############################################

mu <- list()
for(i in 1:length(n_dates))
{
  mu[[i]] <- rep(1.0,n_dates[i]-1)
}
names(mu) <- names(n_dates)
sigma <- mu

theta <- list(zeta = 0.05, # zeta is the probability for a date to be misrecorded, conditional on being recorded (<-> Ei != - 1)
              # TODO:
              # could consider having zeta being type of date specific (e.g. more error on onset than death dates),
              # time specific and/or space specific
              mu = mu, # mean of gamma distributions used to characterise the various delays in different groups: mu[[g]][k] is the mean k^th delay in group g
              sigma = sigma) # sigma of gamma distributions used to characterise the various delays in different groups: sigma[[g]][k] is the std k^th delay in group g


###############################################
### define augmented data ###
###############################################

### D contains the unobserved true dates ###

D <- list() ######### RICH MAYBE IMPROVE THIS BIT OF SHIT CODE ########
for(i in 1:length(n_dates))
{
  D[[i]] <- dat_by_group[[i]]
  for(e in 1:nrow(D[[i]]))
  {
    if(any(is.na(D[[i]][e,])))
    {
      tmp <- which(is.na(D[[i]][e,]))
      if(1 %in% tmp) # dealing with missing values ahead of the series of dates
      {
        min_non_NA_value <- min(which(!is.na(D[[i]][e,])))-1
        for(f in min_non_NA_value:1)
        {
          D[[i]][e,f] <- D[[i]][e,f+1]
        }
      }
      if(any(is.na(D[[i]][e,]))) # dealing with remaining missing values if any
      {
        tmp <- which(is.na(D[[i]][e,]))
        for(f in tmp)
        {
          D[[i]][e,f] <- D[[i]][e,f-1]
        }
      }
    }
  }
}

### E contains an indicator of whether the observed date is the true one or not: ###
### E = -1 if date is unobserved i.e. dat_by_group = -1 ###
### E = 0 if date is observed and exact i.e. dat_by_group = D ###
### E = 1 if date is observed and unexact i.e. dat_by_group not necessarily = D ###

E <- list()
for(i in 1:length(n_dates))
{
  E[[i]] <- as.data.frame(matrix(NA,nrow(dat_by_group[[i]]),ncol(dat_by_group[[i]])))
  for(j in 1:ncol(dat_by_group[[i]]))
  {
    E[[i]][,j] <- rbinom(nrow(dat_by_group[[i]]), 1, theta$zeta)
    E[[i]][is.na(dat_by_group[[i]][,j]),j] <- -1
  }
  names(E[[i]]) <- names(dat_by_group[[i]])
}

### Note that at the moment E can be 1 i.e. there is an error, yet D=dat_by_group i.e. in effect the error is null
### TO DO: TRY AND CHANGE SO THAT E=0 <-> dat_by_group=D

aug_dat <- list(D = D,
                E = E)


###############################################
### likelihood function ###
###############################################

LL_observation_term<-function()
{
  return(1)
}

LL_error_term<-function(E,zeta,log=T)
{
  number_of_errors<-0
  number_of_dates<-0
  for(i in 1:length(E))
  {
    number_of_errors<-number_of_errors+sum(E[[i]]==1)
    number_of_dates<-number_of_dates+prod(dim(E[[i]]))
  }
  result<-dbinom(number_of_errors,number_of_dates,zeta,log=log)

  return(result)
}

LL_delays_term<-function()
{
  return(1)
}

