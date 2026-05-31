#####################################
# Name of the file: dynamic_integrator.R
#####################################
# Description:
# This function models the dynamics of a system with species 
# and resource abundances over time using a set of differential 
# equations.
#####################################
## Author: Karim Ahmed Belmonte
## Email:  kabmisb@gmail.com
## Date:   2024-07-26
#####################################
## INPUT:
## - t: A numeric value representing time.
## - abundances: A numeric vector containing species (S) and resource (R) abundances.
## - parameters: A list
##
## OUTPUT:
## - A list containing a numeric vector of the time derivatives 
##   of species and resource abundances (dSdt and dRdt).
##
## DEPENDENCIES:
## - Requires base R functions (no external packages needed).
##
## USAGE:
## Example:
## params <- list(
##   Cmat = matrix(runif(9), 3, 3),
##   Smat = matrix(runif(9), 3, 3),
##   l = 0.5,
##   w = runif(3),
##   chi = runif(3),
##   g = 0.1,
##   mu = 0.01,
##   tau = 1,
##   kappa = 0.5
## )
## abundances <- runif(6)
## result <- CRm_Modsecret(0, abundances, params)
#####################################

CRm_Modsecret <- function(t, abundances, parameters) {
  
  n_spec    <- parameters$n_spec
  n_res     <- parameters$n_res
  threshold <- parameters$threshold
  abundances[1:n_spec] <- pmin(abundances[1:n_spec], threshold)
  
  S <- abundances[1:n_spec] # Vector that saves abundances of species
  R <- abundances[(n_spec + 1):(n_spec + n_res)] # Vector that saves abundances of resources
  Cmat <- parameters$Cmat
  Smat <- parameters$Smat
  l <- parameters$l
  w <- parameters$w
  chi <- parameters$chi
  g <- parameters$g
  mu <- parameters$mu
  tau <- parameters$tau
  kappa <- parameters$kappa
  
  
  S[S < 1e-8] <- 0 # To prevent numerical problems
  R[R < 1e-8] <- 0 
  
  # Equations
  
  dSdt <- g * (1 - l) * S * (Cmat %*% (R * w) - mu)
  dRdt <- kappa - R / tau - R * (t(Cmat) %*% S) + diag(chi) %*% t(Smat) %*% diag(l * S, n_spec) %*% Cmat %*% (R * w)
  
  return(list(c(dSdt, dRdt)))
}



