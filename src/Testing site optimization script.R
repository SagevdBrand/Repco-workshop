####################################################
####################################################
####### TESTING OPTIMIZATION REGRESSION COEFFICIENTS
set.seed(123)

# Load necessary libraries
source("./config/libraries.R")

# Load functions
source("./src/Functions.R")

# Load data
dm <- readRDS("./data/raw/dm_development.RDS")
dm_val <- readRDS("./data/raw/dm_validation.RDS")

###############
## Define R2 ##
###############

# These values were determined in a different part of my study.
# They were numerically approximated given an AUC of 0.75 and
# event-rate (either 0.05, 0.2 or 0.5).
R2 <- c(0.04131983, 0.12658143, 0.18407039)

#################################
## Define number of predictors ##
#################################
n_pred <- 10

############# OPTIMIZATION PROCESS ############

# Optimizing the dist_R2_prev to determine the beta coefficients
# for which the squared distance between the observed and
# preferred values of R2 and prevalence are minimal.
# par_i are just initial values for the coefficients.

## Through trial and error, these were found to give the best results:
## Initial coefficients for event rate at 0.05, AUC = 0.75 (n_pred = 10)
## -3.37, 0.13
## Initial coefficients for event rate at 0.2, AUC  = 0.75 (n_pred = 10)
## -1.66, 0.13
## Initial coefficients for event rate at 0.5, AUC  = 0.75 (n_pred = 10)
## 0.01, 0.13

##################
### PREV = .05 ###
##################

# Define initial parameters
par1_i <- c(-3.37, 0.13)

# Run optimization script n times, to get different
# values for regression parameters
system.time(results1 <- replicate(n = 5,
                                 optim(par1_i,
                                       dist_R2_prev,
                                       pref_R2 = R2[1],
                                       pref_prev = 0.05,
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence,
# when using the combination of params for each repetition of the above.
results_check1 <- apply(sapply(results1, "[[", 1), 2, checking)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check1, 1, summary)

# If it looks good, take the median of the params from results1
par1 <- apply(sapply(results1, "[[", 1), 1, median)

# Double check whether these values really come down to the
# preferred AUC value and prevalence.
# AUC should be about 0.75 and prevalence 0.05
checking_val(par = par1)

##################
### PREV = 0.2 ###
##################

# Define initial parameters
# These have been found by trial and error
par2_i <- c(-1.66, 0.13)

# Run optimization script n times,
# to get different values for regression parameters
system.time(results2 <- replicate(n = 5,
                                 optim(par2_i,
                                       dist_R2_prev,
                                       pref_R2 = R2[2],
                                       pref_prev = 0.2,
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence,
# when using the combination of params for each repetition of the above.
results_check2 <- apply(sapply(results2, "[[", 1), 2, checking)

# Give a summary of the AUC values (Cstat) and prevalence (prev).
apply(results_check2, 1, summary)

# If it looks good, take the median of the params from results2
par2 <- apply(sapply(results2, "[[", 1), 1, median)

# Double check whether these values really
# come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.2
checking_val(par = par2)

##################
### PREV = 0.5 ###
##################

# Define initial parameters
# These have been found by trial and error
par3_i <- c(0.01, 0.13)

# Run optimization script n times
# to get different values for regression parameters
system.time(results3 <- replicate(n = 5,
                                 optim(par3_i,
                                       dist_R2_prev,
                                       pref_R2 = R2[3],
                                       pref_prev = 0.5,
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence,
# when using the combination of params for each repetition of the above.
results_check3 <- apply(sapply(results3, "[[", 1), 2, checking)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check3, 1, summary)

# If it looks good, take the median of the params from results3
par3 <- apply(sapply(results3, "[[", 1), 1, median)

# Double check whether these values really come
# down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.5
checking_val(par = par3)

#########################
## Saving param values ##
#########################

saveRDS(object = par1, file = "./data/processed/regression_par_prev_0.05.RDS")
saveRDS(object = par2, file = "./data/processed/regression_par_prev_0.2.RDS")
saveRDS(object = par3, file = "./data/processed/regression_par_prev_0.5.RDS")

###################
###################
######## END SCRIPT