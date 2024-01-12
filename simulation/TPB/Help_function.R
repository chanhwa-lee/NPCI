.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(statmod)
library(lme4)
library(SuperLearner)

###---------------------- Help Functions ----------------------###


###-------------- Simulation Data set generation ---------------###
#' Simulation data generation function
#'
#' @param n An integer. Number of clusters
#' @param H A function. Information about treatment model (f_X and sigma)
#' @param G A function. Information about outcome model
#' @return Simulated data following H and G models
#' @examples 
#' data.par = data.sim(1000, H.par.true, G.par.true)
#' data.non = data.sim(1000, H.nonpar.true, G.nonpar.true)

data.sim <- function(n, H, G){
  
  data.list <- list()
  
  N = sample(x = 5:20, size = n, replace = T)
  
  for(i in 1:n){
    
    # Step1. Data generation
    X1 = rnorm(N[i], 0, 1)
    X2 = rbinom(N[i], 1, 0.5)
    C  = rnorm(1, 0, 1)      # Cluster-level covariate
    
    data = data.frame(X1 = X1, X2 = X2, C = C)
    
    # Step2. Treatment model
    p.A = H(data)
    A   = rbinom(N[i], 1, p.A)
    g.A = (sum(A)-A)/(N[i]-1)
    data$A = A
    data$g.A = g.A
    
    # Step3. Outcome model
    p.Y = G(data)
    Y   = rbinom(N[i], 1, p.Y)
    
    data.list[[i]] <- data.frame(id = i,
                                 Y = Y,
                                 A = A,
                                 X1 = X1,
                                 X2 = X2,
                                 C = C)
  }
  
  return(dplyr::bind_rows(data.list))
}


###       Parametric DGP         ###
G.par.true = function(data){
  return(with(data, plogis(3 - 2*A - g.A - 1.5*X1 + 2*X2 - 2*C)))
}

H.par.true = function(data){
  return(with(data, plogis(0.1 + 0.2*X1 + 0.2*X2 + 0.1*C)))
}

###      NonParametric DGP      ###
G.nonpar.true = function(data){
  return(with(data, plogis(3 - 2*A - g.A - 1.5*abs(X1) + 2*X2 - 3*abs(X1)*X2 - 2*I(C>0))))
}

H.nonpar.true = function(data){
  return(with(data, plogis(0.1 + 0.2*abs(X1) + 0.2*abs(X1)*X2 + 0.1*I(C>0))))
}