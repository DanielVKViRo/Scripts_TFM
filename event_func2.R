#####################################
# Name of the file: convergence_param_func.R
#####################################
# Description: Script for a function that will calculate the parameters of the 
# convergence of the system
#####################################
## Author: Karim Ahmed Belmonte
## Email:  kabmisb@gmail.com
## Date:     2024-09-24
#####################################
## INPUT: Abundinit, times
## OUTPUT:
## DEPENDENCIES:
## USAGE: 
#####################################

event_func2 <- function(t, y, parms) {
  #browser()
  # Define an event function to stop integration when convergence is achieved
  t_inside <<- t_inside + 1 # Marks the integration step
  convergence_condition_met <- 0
  new_abund <- y
  if(t_inside > n_steps + 1){
    #browser()
    y_mean <- colMeans(s_r_global[1:(n_steps-1),])
    y_star <- new_abund
    #y_init <- s_r_global[,1]
    epsilon <- abs(y_mean - y_star) #/y_init # could be infty for initially empty compartments
    if (all(epsilon < thr)) {
      key_conv <<- key_conv + 1
    } else {
      key_conv <<- 0
    }
  }
  s_r_index <<- s_r_index + 1
  if ((t_inside %% n_steps) == 0 ){
    s_r_index <<- 0
  }
  s_r_global[s_r_index, ] <<- y
  
  if (key_conv == n_conv + 1 ){ # convergence condition was achieved n_conv times
    convergence_condition_met = 1
  }
  # Check for convergence condition
  if (convergence_condition_met == 1) {
    return(0)  # Return 0 to stop integration
  } else {
    return(-1)   # Return -1 to continue integration
  }
}