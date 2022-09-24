# load necessary packages
library(tidyverse)
library(mvtnorm)
library(mgcv)
library(deSolve)

Probabilities_all_Alf <- function(alpha,S,sims,delta_t) {
  survival <- rep(0,S)
  # sims <- 250  ### number of simualtions to obtain the probabilities of species persistence. In the text we used 100,000
  y <- 1
  while(y <= sims){
    time_step <- seq(0,300,by=delta_t)
    N0 <- rep(1,S) #initial conditions for species abundances
    N0 <- N0 / sum(N0)
    r <- sphere_sampling_ALF(S) ### sampling random K's in the unit sphere
    parms <- list(r=r, alpha = alpha) ##ODE
    model <- function(t,N,parms){ dN <- N * (parms$r - parms$alpha %*% N); list(dN)}
    sol <- ode(N0,time_step,model,parms,
               method = 'bdf', 
               rtol = 1e-15, 
               atol = 1e-15, 
               maxsteps = 1000000)
    
    #plot(sol)
    
    survival_aux <- rep(0,S)
    
    for(z in 1:S){
      survival_aux[z] <- (sol[nrow(sol),z+1]>0.000001)*1 ### check if the species survived in the simulation
    }
    
    if(!all(is.na(survival_aux))){
      survival <- survival + survival_aux
      y <- y + 1
    }
    cat(y,"\n")
    
  }
  survival <- survival / sims ### calculate probabilities
  out <- survival
  return(out)
}

Probabilities_all_Alf_POSITIVE <- function(alpha,S,sims,delta_t) {
  survival <- rep(0,S)
  # sims <- 250  ### number of simualtions to obtain the probabilities of species persistence. In the text we used 100,000
  y <- 1
  while(y <= sims){
    time_step <- seq(0,300,by=delta_t)
    N0 <- rep(1,S) #initial conditions for species abundances
    N0 <- N0 / sum(N0)
    r <- abs(sphere_sampling_ALF(S)) ### sampling random K's in the unit sphere
    parms <- list(r=r, alpha = alpha) ##ODE
    model <- function(t,N,parms){ dN <- N * (parms$r - parms$alpha %*% N); list(dN)}
    sol <- ode(N0,time_step,model,parms)
    
    #plot(sol)
    
    survival_aux <- rep(0,S)
    
    for(z in 1:S){
      survival_aux[z] <- (sol[nrow(sol),z+1]>0.000001)*1 ### check if the species survived in the simulation
    }
    
    if(!all(is.na(survival_aux))){
      survival <- survival + survival_aux
      y <- y + 1
    }
    cat(y,"\n")
    
  }
  survival <- survival / sims ### calculate probabilities
  out <- survival
  return(out)
}

# Samples m vectors randomly on a n-dimensional unit sphere
sphere_sampling_ALF <- function(m) {
  r <- (rnorm(m))
  d <- sqrt(sum(r^2))
  r <- r/d
  return(r)
}