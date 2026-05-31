#####################################
# Name of the file: convergence_param_func2.R
#####################################
# Description: 
# Script for a function that calculates the parameters of the 
# convergence of the system.
#####################################
## Author: Karim Ahmed Belmonte
## Email:  kabmisb@gmail.com
## Date:   2024-10-04
#####################################
## INPUT:
## - Abundinit: A numeric vector representing initial species abundances.
## - times: A numeric value representing the time steps for the convergence check.
##
## OUTPUT:
## - The function initializes and sets global parameters for convergence checks.
##   However, it currently does not return a value.
##
## DEPENDENCIES:
## - Base R functions (no external packages required).
##
## USAGE:
## Example:
## Abundinit <- runif(10, 0, 1) # Random initial abundances
## times <- 100 # Number of time steps
## convergence_param_func2(Abundinit, times)
#####################################

convergence_param_func2 <-  function(Abundinit, times) {
  
  #browser() #debug
  TOL <<- 1e-2
  thr <<- rep(TOL, length(Abundinit)) # threshold resources, species, and convergence
  t_inside <<- 0
  n_steps <<- 51 # minimum number of steps to determine convergence
  s_r_global <<- matrix(NA, nrow = n_steps, ncol = length(Abundinit))
  s_r_index <<- 0
 
  n_conv <<- 25 # number of times convergence should be verified to halt
  key_conv <<- 0
  
  return}
