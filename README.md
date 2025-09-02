# MixDiff

MixDiff is an R package for Bayesian estimation of event dates and delay
distributions using epidemiological line-list data. It is designed to handle
commonly encountered issues with line-list datasets during outbreaks, including
missing event dates (e.g., symptom onset, hospitalisation, or death) and errors
in recorded dates.

To do:
- Finish tests and set up automated testing
- Break code down into functions to test more easily
- Add more checks for inputs
- Document functions
- Add vignettes and link to them here
- Decide on package name
- Decide on other names e.g. index_dates --> delays, move_D --> move_dates,
    move_E --> move_errors
- Decide on input format: single input dataset with dates and group variables
    --> user would need to distinguish between missing vs NA depending on group
- Change index_dates/delays format to take date names
    --> table with "from", "to" and "group" columns
    --> include hyperparameters in this table?
- Change fraction_move_Di?
- Function to set up MCMC_settings object (like make_config()?)
- How to handle delays shared between groups - single delay/estimate per group/both