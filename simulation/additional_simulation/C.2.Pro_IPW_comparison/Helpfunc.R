#*******************************************************************************
#**********           C.2. Proposed vs IPW comparison                 **********
#**********           under CIPS policy                              **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********                                                           **********
#**********           Version: 2.0                                    **********
#**********           Jan 12, 2024                                    **********
#*******************************************************************************

### Libraries ###
library(dplyr)


###----------- `deltas` of interest ---------------###

deltas = c(0.5,1,2)


###----------- Cluster size (N) distribution ---------------###

N.dist = function(m) sample(x = 5:20, size = m, replace = T)


###----------- Estimand computation function ---------------###

#' Causal estimands computation under the Incremental Propensity Score policy
#'
#' @param N An integer. Cluster size
#' @param deltas A numeric vector. Deltas of IPS policies. Should include 1 
#' (delta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Approximated causal estimands under the IPS policy over the delta values

DGP <- function(N, type = "data", deltas = NULL, r = NULL){
  
  ### Observed data generation ###

  # Step1. Covariate generation
  X1 = rnorm(N, 0, 1)
  X2 = rbinom(N, 1, 0.5)
  C  = rnorm(1, 0, 1)

  # Step2. Treatment model
  pi <- plogis(0.1 + 0.2*abs(X1) + 0.2*abs(X1)*X2 + 0.1*I(C>0))
  A = rbinom(N, 1, pi)
  g.A = (sum(A)-A)/(N-1)

  # Step3. Outcome model
  p.Y = function() plogis(3 - 2*A - g.A - 1.5*abs(X1) + 2*X2 - 3*abs(X1)*X2 - 2*I(C>0))
  Y = rbinom(N, 1, p.Y())

  if(type == "data"){
    
    ### Observed data ###
    
    data <- data.frame(Y = Y, A = A, 
                       X1 = X1, X2 = X2, C = C)
    return(data)
  }
  
  if(type == "estimands"){
      
    #### Estimand computation over `delta` values ###
    
    # Add `delta` = 1 if not in `deltas`
    if(!(1 %in% deltas)) deltas <- c(deltas, 1)
    
    estimands <- data.frame(delta = deltas,
                            mu = 0, mu_1 = 0, mu_0 = 0)
    
    for(i in seq_len(length(deltas))){
      
      delta = deltas[i]
      
      pi.delta <- delta*pi / (delta*pi + 1 - pi)
      
      ### Subsampling approximation ###
      A <- rbinom(n = r*N, size = 1, prob = rep(pi.delta, r))
      
      g.A = (rep(aggregate(A, list(rep(1:r, each = N)), sum)[,-1], each = N) - A) / (N-1)
      
      mu   = mean(p.Y(), na.rm = T)
      
      mu_1 = mean(p.Y() * A / rep(pi.delta, r), na.rm = T)
      
      mu_0 = mean(p.Y() * (1-A) / rep(1-pi.delta, r), na.rm = T)
      
      estimands[i, ] = c(delta, mu, mu_1, mu_0)
      
    }
    
    return(estimands)
  
  }
  
  if(type == "alpha.comp"){
    
    ### `alpha` corresponding to `delta` for Barkley method ###
    
    alphas <- data.frame(delta = deltas,
                         alpha = 0)
    
    for(i in seq_len(length(deltas))){
      
      delta = deltas[i]
      
      pi.delta <- delta*pi / (delta*pi + 1 - pi)
      
      alphas[i, "alpha"] = mean(pi.delta)
      
    }
    
    return(alphas)
    
  }
  
}



###----------- Observed data simulation function ---------------###

#' Causal estimands computation under the Incremental Propensity Score policy
#'
#' @param N An integer. Cluster size
#' @param deltas A numeric vector. Deltas of IPS policies. Should include 1 
#' (delta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Approximated causal estimands under the IPS policy over the delta values

data.sim = function(m){
  N.list = N.dist(m)
  data = lapply(N.list, DGP)
  data = dplyr::bind_rows(data,.id = "id") %>% mutate(id = as.numeric(id))
  return(data)
}

###----------- Estimands simulation function ---------------###

#' Causal estimands computation under the Incremental Propensity Score policy
#'
#' @param N An integer. Cluster size
#' @param deltas A numeric vector. Deltas of IPS policies. Should include 1 
#' (delta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Approximated causal estimands under the IPS policy over the delta values

estimands.sim = function(m, deltas, r){
  N = N.dist(m)
  estimands = lapply(N, DGP, type = "estimands", deltas = deltas, r = r)
  estimands = dplyr::bind_rows(estimands) %>% 
    group_by(delta) %>% 
    summarise(mu = mean(mu),
              mu_1 = mean(mu_1),
              mu_0 = mean(mu_0))
  
  ### Causal effects computation using delta==1 as the standard ###
  standard = estimands %>% filter(delta == 1)
  
  estimands = estimands %>% mutate(de  = mu_1 - mu_0, 
                                   se_1 = mu_1 - standard$mu_1,
                                   se_0 = mu_0 - standard$mu_0,
                                   oe  = mu   - standard$mu,
                                   te  = mu_1 - standard$mu_0)
  
  return(estimands)
}


###-------------- Estimator (Proposed vs IPW) function --------------------###

estimator.simul <- function(dat, X.trt, X.out, deltas, nsplits, r){
  
  ## If deltas does not include 1, then add it to the list
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  
  ## Split dat into lists to get N for each cluster
  dat.list <- split(dat, f = dat$id)
  
  ## setup storage
  m <- length(dat.list)        # Number of clusters
  N <- sapply(dat.list, nrow)  # Number of individuals in each cluster
  N.limit <- cumsum(c(0,N))
  
  est.par <- se.par <-  data.frame(delta = deltas,
                                   mu = 0, mu_1 = 0, mu_0 = 0,
                                   de = 0, se_1 = 0, se_0 = 0,
                                   oe = 0, te = 0)
  
  est.nonpar <- se.nonpar <-  data.frame(delta = deltas,
                                         mu = 0, mu_1 = 0, mu_0 = 0,
                                         de = 0, se_1 = 0, se_0 = 0,
                                         oe = 0, te = 0)
  
  est.par.ipw <- se.par.ipw <-  data.frame(delta = deltas,
                                           mu = 0, mu_1 = 0, mu_0 = 0,
                                           de = 0, se_1 = 0, se_0 = 0,
                                           oe = 0, te = 0)
  
  est.nonpar.ipw <- se.nonpar.ipw <-  data.frame(delta = deltas,
                                                 mu = 0, mu_1 = 0, mu_0 = 0,
                                                 de = 0, se_1 = 0, se_0 = 0,
                                                 oe = 0, te = 0)
  
  pi.par <- numeric(sum(N))
  mu.par <- numeric(sum(N))
  
  pi.nonpar <- numeric(sum(N))
  mu.nonpar <- numeric(sum(N))
  
  ifvals = list(mu   = matrix(0, nrow = length(deltas), ncol = m),
                mu_1 = matrix(0, nrow = length(deltas), ncol = m),
                mu_0 = matrix(0, nrow = length(deltas), ncol = m))
  
  rownames(ifvals$mu) = deltas
  rownames(ifvals$mu_1) = deltas
  rownames(ifvals$mu_0) = deltas
  
  ifvals.list = list(pro.p  = ifvals,
                     pro.np = ifvals,
                     ipw.p  = ifvals,
                     ipw.np = ifvals)
  
  fold <- sample(1:m, replace = F) %% nsplits + 1
  foldlong <- rep(fold, times = N)
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:nsplits) {
    
    print(paste("   Split:", split))
    
    train.idx = which(foldlong != split)
    eval.idx  = which(foldlong == split)
    
    ### 1. Main-effect only GLM (Logistic regression)
    pi.fit.par <- glm(A ~ ., data = cbind(A = dat$A, X.trt)[train.idx,], family = "binomial")
    
    pi.par[eval.idx] <- predict(pi.fit.par, newdata = X.trt[eval.idx,], type = "response")
    
    mu.fit.par <- glm(Y ~ ., data = cbind(Y = dat$Y, A = dat$A, g.A = dat$g.A, X.out)[train.idx,], family = "binomial")
    
    mu.par[eval.idx] <- predict(mu.fit.par, newdata = cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], type = "response")
    
    ### 2. SuperLearner (RF & GLM)
    pi.fit.nonpar = SuperLearner(Y = dat$A[train.idx], X = X.trt[train.idx,],
                                 family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    pi.nonpar[eval.idx] <- predict(pi.fit.nonpar, X.trt[eval.idx,], onlySL = TRUE)$pred
    
    mu.fit.nonpar <- SuperLearner(Y = dat$Y[train.idx], X = cbind(A = dat$A, g.A = dat$g.A, X.out)[train.idx,],
                                  family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    mu.nonpar[eval.idx] <- predict(mu.fit.nonpar, cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], onlySL = TRUE)$pred
    
    ## evaluate IF values for each cluster in eval fold
    for(i in which(fold == split)){
      
      # print(paste("   Cluster:", i))
      
      dat.sub.idx <- (N.limit[i]+1):N.limit[i+1]
      dat.sub <- dat[dat.sub.idx,]
      
      pi.sub.par  <- pi.par[dat.sub.idx]
      mu.sub.par  <- mu.par[dat.sub.idx]
      
      pi.sub.nonpar  <- pi.nonpar[dat.sub.idx]
      mu.sub.nonpar  <- mu.nonpar[dat.sub.idx]
      
      # random binary vector from A(N_i) w. uniform sampling, instead of summing all vectors in A(N)
      
      # Help function
      agg.fun <- function(a, fun){
        aggregate(a, list(rep(1:r, each = N[i])), fun)[,-1]
      }
      
      a <- sample(c(0,1), r*N[i], replace = T)                                   # rN[i] vector
      g.a = (rep(agg.fun(a, sum), each = N[i]) - a) / (N[i]-1)
      
      
      ### 1. Parametric
      
      ## Outcome regression
      mu.a.par <- predict(mu.fit.par,
                          newdata = cbind(data.frame(A = a, g.A = g.a),
                                          bind_rows(replicate(r, X.out[dat.sub.idx,], simplify = FALSE))),
                          type = "response")
      
      pi.delta.sub.par <- pi.sub.par %o% deltas / (pi.sub.par %o% deltas + 1 - pi.sub.par)       # N[i] x length(deltas) matrix
      
      cluster.prob.ratio.par <- deltas^(sum(dat.sub$A)) / 
        apply(pi.sub.par %o% deltas + 1 - pi.sub.par, 2, prod)                           # length(deltas) vector
      
      pi.delta.sub.rep.par <- rep(1, r) %x% pi.delta.sub.par                             # rN[i] x length(deltas) matrix (replicate pi.delta.sub r times)
      
      cluster.prob.par <- agg.fun(a*pi.delta.sub.rep.par + (1-a)*(1-pi.delta.sub.rep.par), prod)   # r x length(deltas) matrix
      
      add.factor.par <-
        (2*a-1) %o% deltas * rep(dat.sub$A - pi.sub.par, r) / 
        (rep(1, r) %x% (pi.sub.par %o% deltas + 1 - pi.sub.par)) /
        (a*rep(pi.sub.par, r) %o% deltas + (1-a)*(1-rep(pi.sub.par, r)))                 # rN[i] x length(deltas) matrix
      
      out.reg.mu.par <- apply(agg.fun(mu.a.par, mean) * (1 + agg.fun(add.factor.par, sum)) * 
                                2^(N[i]) * cluster.prob.par, 2, mean)                    # length(deltas) vector
      
      out.reg.mu_1.par <-
        apply(agg.fun( (mu.a.par * a) %o% rep(1, length(deltas)) / pi.delta.sub.rep.par * 
                         (1 + as.matrix(agg.fun(add.factor.par, sum)) %x% rep(1,N[i]) - add.factor.par), mean) *
                2^(N[i]) * cluster.prob.par, 2, mean)                                 # length(deltas) vector
      
      out.reg.mu_0.par <-
        apply(agg.fun( (mu.a.par * (1-a)) %o% rep(1, length(deltas)) / (1-pi.delta.sub.rep.par) * 
                         (1 + as.matrix(agg.fun(add.factor.par, sum)) %x% rep(1,N[i]) - add.factor.par), mean) *
                2^(N[i]) * cluster.prob.par, 2, mean)                                # length(deltas) vector
      
      ## Bias correction part for mu
      bias.cor.mu.par <- mean(dat.sub$Y - mu.sub.par, na.rm = T) * cluster.prob.ratio.par    # length(deltas) vector
      
      ## Bias correction part for mu_1
      bias.cor.mu_1.par <-
        apply( (dat.sub$Y - mu.sub.par) * I(dat.sub$A==1) / pi.delta.sub.par, 2, mean) * 
        cluster.prob.ratio.par                                                       # length(deltas) vector
      
      ## Bias correction part for mu_0
      bias.cor.mu_0.par <-
        apply( (dat.sub$Y - mu.sub.par) * I(dat.sub$A==0) / (1-pi.delta.sub.par), 2, mean) * 
        cluster.prob.ratio.par                                                       # length(deltas) vector
      
      ## IF for mu at cluster i
      ifvals.list$pro.p$mu[,i] <- out.reg.mu.par + bias.cor.mu.par
      ifvals.list$ipw.p$mu[,i] <- mean(dat.sub$Y, na.rm = T) * cluster.prob.ratio.par
      
      ## IF for mu_1 at cluster i
      ifvals.list$pro.p$mu_1[,i] <- out.reg.mu_1.par + bias.cor.mu_1.par
      ifvals.list$ipw.p$mu_1[,i] <- 
        apply(dat.sub$Y * I(dat.sub$A==1) / pi.delta.sub.par, 2, mean) * cluster.prob.ratio.par
      
      ## IF for mu_0 at cluster i
      ifvals.list$pro.p$mu_0[,i] <- out.reg.mu_0.par + bias.cor.mu_0.par
      ifvals.list$ipw.p$mu_0[,i] <-
        apply(dat.sub$Y * I(dat.sub$A==0) / (1-pi.delta.sub.par), 2, mean) * cluster.prob.ratio.par
      
      
      ### 2. Nonparametric
      
      ## Outcome regression
      mu.a.nonpar <- predict(mu.fit.nonpar, 
                             cbind(data.frame(A = a, g.A = g.a), bind_rows(replicate(r, X.out[dat.sub.idx,], simplify = FALSE))),
                             onlySL = TRUE)$pred
      
      
      pi.delta.sub.nonpar <- pi.sub.nonpar %o% deltas / (pi.sub.nonpar %o% deltas + 1 - pi.sub.nonpar)       # N[i] x length(deltas) matrix
      
      cluster.prob.ratio.nonpar <- deltas^(sum(dat.sub$A)) / 
        apply(pi.sub.nonpar %o% deltas + 1 - pi.sub.nonpar, 2, prod)                           # length(deltas) vector
      
      pi.delta.sub.rep.nonpar <- rep(1, r) %x% pi.delta.sub.nonpar                             # rN[i] x length(deltas) matrix (replicate pi.delta.sub r times)
      
      cluster.prob.nonpar <- agg.fun(a*pi.delta.sub.rep.nonpar + (1-a)*(1-pi.delta.sub.rep.nonpar), prod)   # r x length(deltas) matrix
      
      add.factor.nonpar <-
        (2*a-1) %o% deltas * rep(dat.sub$A - pi.sub.nonpar, r) / 
        (rep(1, r) %x% (pi.sub.nonpar %o% deltas + 1 - pi.sub.nonpar)) /
        (a*rep(pi.sub.nonpar, r) %o% deltas + (1-a)*(1-rep(pi.sub.nonpar, r)))                 # rN[i] x length(deltas) matrix
      
      out.reg.mu.nonpar <- apply(agg.fun(mu.a.nonpar, mean) * (1 + agg.fun(add.factor.nonpar, sum)) * 
                                   2^(N[i]) * cluster.prob.nonpar, 2, mean)                    # length(deltas) vector
      
      out.reg.mu_1.nonpar <-
        apply(agg.fun( (mu.a.nonpar * a)[,rep(1, length(deltas))] / pi.delta.sub.rep.nonpar * 
                         (1 + as.matrix(agg.fun(add.factor.nonpar, sum)) %x% rep(1,N[i]) - add.factor.nonpar), mean) *
                2^(N[i]) * cluster.prob.nonpar, 2, mean)                                 # length(deltas) vector
      
      out.reg.mu_0.nonpar <-
        apply(agg.fun( (mu.a.nonpar * (1-a))[,rep(1, length(deltas))] / (1-pi.delta.sub.rep.nonpar) * 
                         (1 + as.matrix(agg.fun(add.factor.nonpar, sum)) %x% rep(1,N[i]) - add.factor.nonpar), mean) *
                2^(N[i]) * cluster.prob.nonpar, 2, mean)                                # length(deltas) vector
      
      ## Bias correction part for mu
      bias.cor.mu.nonpar <- mean(dat.sub$Y - mu.sub.nonpar, na.rm = T) * cluster.prob.ratio.nonpar    # length(deltas) vector
      
      ## Bias correction part for mu_1
      bias.cor.mu_1.nonpar <-
        apply( (dat.sub$Y - mu.sub.nonpar) * I(dat.sub$A==1) / pi.delta.sub.nonpar, 2, mean) * 
        cluster.prob.ratio.nonpar                                                       # length(deltas) vector
      
      ## Bias correction part for mu_0
      bias.cor.mu_0.nonpar <-
        apply( (dat.sub$Y - mu.sub.nonpar) * I(dat.sub$A==0) / (1-pi.delta.sub.nonpar), 2, mean) * 
        cluster.prob.ratio.nonpar                                                       # length(deltas) vector
      
      
      
      ## IF for mu at cluster i
      ifvals.list$pro.np$mu[,i] <- out.reg.mu.nonpar + bias.cor.mu.nonpar
      ifvals.list$ipw.np$mu[,i] <- mean(dat.sub$Y, na.rm = T) * cluster.prob.ratio.nonpar
      
      ## IF for mu_1 at cluster i
      ifvals.list$pro.np$mu_1[,i] <- out.reg.mu_1.nonpar + bias.cor.mu_1.nonpar
      ifvals.list$ipw.np$mu_1[,i] <- 
        apply(dat.sub$Y * I(dat.sub$A==1) / pi.delta.sub.nonpar, 2, mean) * cluster.prob.ratio.nonpar
      
      ## IF for mu_0 at cluster i
      ifvals.list$pro.np$mu_0[,i] <- out.reg.mu_0.nonpar + bias.cor.mu_0.nonpar
      ifvals.list$ipw.np$mu_0[,i] <-
        apply(dat.sub$Y * I(dat.sub$A==0) / (1-pi.delta.sub.nonpar), 2, mean) * cluster.prob.ratio.nonpar
      
    }
    
  }
  
  ## helper function
  IF.to.est <- function(IF, fold){
    est = mean(aggregate(IF, by = list(fold), mean)[,-1])
    se  = sqrt(1/m * (mean(aggregate(IF, by = list(fold), function(phi) mean(phi^2))[,-1]) - est^2))
    return(c(est, se))
  }
  
  
  IF.list.to.est <- function(IF.list, fold){
    
    est <- se <-  data.frame(delta = deltas,
                             mu = 0, mu_1 = 0, mu_0 = 0,
                             de = 0, se_1 = 0, se_0 = 0,
                             oe = 0, te = 0)
    
    for(delta.idx in seq_len(length(deltas))){
      
      # delta specific estimates and se estimates
      est[delta.idx, "mu"]   <- IF.to.est(IF.list$mu  [delta.idx,], fold)[1]
      est[delta.idx, "mu_1"] <- IF.to.est(IF.list$mu_1[delta.idx,], fold)[1]
      est[delta.idx, "mu_0"] <- IF.to.est(IF.list$mu_0[delta.idx,], fold)[1]
      est[delta.idx, "de"]   <- IF.to.est(IF.list$mu_1[delta.idx,] - IF.list$mu_0[delta.idx,], fold)[1]
      
      se[delta.idx, "mu"]   <- IF.to.est(IF.list$mu  [delta.idx,], fold)[2]
      se[delta.idx, "mu_1"] <- IF.to.est(IF.list$mu_1[delta.idx,], fold)[2]
      se[delta.idx, "mu_0"] <- IF.to.est(IF.list$mu_0[delta.idx,], fold)[2]
      se[delta.idx, "de"]   <- IF.to.est(IF.list$mu_1[delta.idx,] - IF.list$mu_0[delta.idx,], fold)[2]
      
      # delta versus delta0 == 1 comparison estimates and se estimates
      delta0.idx = which(deltas == 1)
      est[delta.idx, "se_1"] <- IF.to.est(IF.list$mu_1[delta.idx,] - IF.list$mu_1[delta0.idx,], fold)[1]
      est[delta.idx, "se_0"] <- IF.to.est(IF.list$mu_0[delta.idx,] - IF.list$mu_0[delta0.idx,], fold)[1]
      est[delta.idx, "oe"]   <- IF.to.est(IF.list$mu  [delta.idx,] - IF.list$mu  [delta0.idx,], fold)[1]
      est[delta.idx, "te"]   <- IF.to.est(IF.list$mu_1[delta.idx,] - IF.list$mu_0[delta0.idx,], fold)[1]
      
      se[delta.idx, "se_1"] <- IF.to.est(IF.list$mu_1[delta.idx,] - IF.list$mu_1[delta0.idx,], fold)[2]
      se[delta.idx, "se_0"] <- IF.to.est(IF.list$mu_0[delta.idx,] - IF.list$mu_0[delta0.idx,], fold)[2]
      se[delta.idx, "oe"]   <- IF.to.est(IF.list$mu  [delta.idx,] - IF.list$mu  [delta0.idx,], fold)[2]
      se[delta.idx, "te"]   <- IF.to.est(IF.list$mu_1[delta.idx,] - IF.list$mu_0[delta0.idx,], fold)[2]
      
    }
    
    return(list(est = est, se = se))
  }
  
  methods = names(ifvals.list)
  estimate =
    lapply(methods,
           function(method) {
             if (!(method %in% c("pro.p", "pro.np"))) fold = rep(1, m)
             return(IF.list.to.est(ifvals.list[[method]], fold))
           })
  names(estimate) = methods
  
  return(estimate)
  
}



###-------------- Proposed Estimator function --------------------###

estimator <- function(dat, X.trt, X.out, deltas, nsplits, r){
  
  ## If deltas does not include 1, then add it to the list
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  
  ## Split dat into lists to get N for each cluster
  dat.list <- split(dat, f = dat$id)
  
  ## setup storage
  m <- length(dat.list)        # Number of clusters
  N <- sapply(dat.list, nrow)  # Number of individuals in each cluster
  N.limit <- cumsum(c(0,N))
  
  est <- se <-  data.frame(delta = deltas,
                                         mu = 0, mu_1 = 0, mu_0 = 0,
                                         de = 0, se_1 = 0, se_0 = 0,
                                         oe = 0, te = 0)
  
  pi <- numeric(sum(N))
  mu <- numeric(sum(N))
  
  ifvals = list(mu   = matrix(0, nrow = length(deltas), ncol = m),
                mu_1 = matrix(0, nrow = length(deltas), ncol = m),
                mu_0 = matrix(0, nrow = length(deltas), ncol = m))
  
  rownames(ifvals$mu) = deltas
  rownames(ifvals$mu_1) = deltas
  rownames(ifvals$mu_0) = deltas
  
  fold <- sample(1:m, replace = F) %% nsplits + 1
  foldlong <- rep(fold, times = N)
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:nsplits) {
    
    print(paste("   Split:", split))
    
    train.idx = which(foldlong != split)
    eval.idx  = which(foldlong == split)
    
    ### SuperLearner (RF & GLM) ###
    pi.fit = SuperLearner(Y = dat$A[train.idx], X = X.trt[train.idx,],
                                 family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    pi[eval.idx] <- predict(pi.fit, X.trt[eval.idx,], onlySL = TRUE)$pred
    
    mu.fit <- SuperLearner(Y = dat$Y[train.idx], X = cbind(A = dat$A, g.A = dat$g.A, X.out)[train.idx,],
                                  family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    mu[eval.idx] <- predict(mu.fit, cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], onlySL = TRUE)$pred
    
    ## evaluate IF values for each cluster in eval fold
    for(i in which(fold == split)){
      
      # print(paste("   Cluster:", i))
      
      dat.sub.idx <- (N.limit[i]+1):N.limit[i+1]
      dat.sub <- dat[dat.sub.idx,]
      
      pi.sub  <- pi[dat.sub.idx]
      mu.sub  <- mu[dat.sub.idx]
      
      # random binary vector from A(N_i) w. uniform sampling, instead of summing all vectors in A(N)
      
      # Help function
      agg.fun <- function(a, fun){
        aggregate(a, list(rep(1:r, each = N[i])), fun)[,-1]
      }
      
      a <- sample(c(0,1), r*N[i], replace = T)                                   # rN[i] vector
      g.a = (rep(agg.fun(a, sum), each = N[i]) - a) / (N[i]-1)
      
      ## Outcome regression
      mu.a <- predict(mu.fit, 
                             cbind(data.frame(A = a, g.A = g.a), bind_rows(replicate(r, X.out[dat.sub.idx,], simplify = FALSE))),
                             onlySL = TRUE)$pred
      
      
      pi.delta.sub <- pi.sub %o% deltas / (pi.sub %o% deltas + 1 - pi.sub)       # N[i] x length(deltas) matrix
      
      cluster.prob.ratio <- deltas^(sum(dat.sub$A)) / 
        apply(pi.sub %o% deltas + 1 - pi.sub, 2, prod)                           # length(deltas) vector
      
      pi.delta.sub.rep <- rep(1, r) %x% pi.delta.sub                             # rN[i] x length(deltas) matrix (replicate pi.delta.sub r times)
      
      cluster.prob <- agg.fun(a*pi.delta.sub.rep + (1-a)*(1-pi.delta.sub.rep), prod)   # r x length(deltas) matrix
      
      add.factor <-
        (2*a-1) %o% deltas * rep(dat.sub$A - pi.sub, r) / 
        (rep(1, r) %x% (pi.sub %o% deltas + 1 - pi.sub)) /
        (a*rep(pi.sub, r) %o% deltas + (1-a)*(1-rep(pi.sub, r)))                 # rN[i] x length(deltas) matrix
      
      out.reg.mu <- apply(agg.fun(mu.a, mean) * (1 + agg.fun(add.factor, sum)) * 
                                   2^(N[i]) * cluster.prob, 2, mean)                    # length(deltas) vector
      
      out.reg.mu_1 <-
        apply(agg.fun( (mu.a * a)[,rep(1, length(deltas))] / pi.delta.sub.rep * 
                         (1 + as.matrix(agg.fun(add.factor, sum)) %x% rep(1,N[i]) - add.factor), mean) *
                2^(N[i]) * cluster.prob, 2, mean)                                 # length(deltas) vector
      
      out.reg.mu_0 <-
        apply(agg.fun( (mu.a * (1-a))[,rep(1, length(deltas))] / (1-pi.delta.sub.rep) * 
                         (1 + as.matrix(agg.fun(add.factor, sum)) %x% rep(1,N[i]) - add.factor), mean) *
                2^(N[i]) * cluster.prob, 2, mean)                                # length(deltas) vector
      
      ## Bias correction part for mu
      bias.cor.mu <- mean(dat.sub$Y - mu.sub, na.rm = T) * cluster.prob.ratio    # length(deltas) vector
      
      ## Bias correction part for mu_1
      bias.cor.mu_1 <-
        apply( (dat.sub$Y - mu.sub) * I(dat.sub$A==1) / pi.delta.sub, 2, mean) * 
        cluster.prob.ratio                                                       # length(deltas) vector
      
      ## Bias correction part for mu_0
      bias.cor.mu_0 <-
        apply( (dat.sub$Y - mu.sub) * I(dat.sub$A==0) / (1-pi.delta.sub), 2, mean) * 
        cluster.prob.ratio                                                       # length(deltas) vector
      
      
      
      ## IF for mu at cluster i
      ifvals$mu[,i] <- out.reg.mu + bias.cor.mu
      
      ## IF for mu_1 at cluster i
      ifvals$mu_1[,i] <- out.reg.mu_1 + bias.cor.mu_1
      
      ## IF for mu_0 at cluster i
      ifvals$mu_0[,i] <- out.reg.mu_0 + bias.cor.mu_0
      
    }
    
  }
  
  ## helper function
  IF.to.est <- function(IF, fold){
    est = mean(aggregate(IF, by = list(fold), mean)[,-1])
    se  = sqrt(1/m * (mean(aggregate(IF, by = list(fold), function(phi) mean(phi^2))[,-1]) - est^2))
    return(c(est, se))
  }
  
  ## Compute estimate and se
  
  for(delta.idx in seq_len(length(deltas))){
    
    # delta specific estimates and se estimates
    est[delta.idx, "mu"]   <- IF.to.est(ifvals$mu  [delta.idx,], fold)[1]
    est[delta.idx, "mu_1"] <- IF.to.est(ifvals$mu_1[delta.idx,], fold)[1]
    est[delta.idx, "mu_0"] <- IF.to.est(ifvals$mu_0[delta.idx,], fold)[1]
    est[delta.idx, "de"]   <- IF.to.est(ifvals$mu_1[delta.idx,] - ifvals$mu_0[delta.idx,], fold)[1]
    
    se[delta.idx, "mu"]   <- IF.to.est(ifvals$mu  [delta.idx,], fold)[2]
    se[delta.idx, "mu_1"] <- IF.to.est(ifvals$mu_1[delta.idx,], fold)[2]
    se[delta.idx, "mu_0"] <- IF.to.est(ifvals$mu_0[delta.idx,], fold)[2]
    se[delta.idx, "de"]   <- IF.to.est(ifvals$mu_1[delta.idx,] - ifvals$mu_0[delta.idx,], fold)[2]
    
    # delta versus delta0 == 1 comparison estimates and se estimates
    delta0.idx = which(deltas == 1)
    est[delta.idx, "se_1"] <- IF.to.est(ifvals$mu_1[delta.idx,] - ifvals$mu_1[delta0.idx,], fold)[1]
    est[delta.idx, "se_0"] <- IF.to.est(ifvals$mu_0[delta.idx,] - ifvals$mu_0[delta0.idx,], fold)[1]
    est[delta.idx, "oe"]   <- IF.to.est(ifvals$mu  [delta.idx,] - ifvals$mu  [delta0.idx,], fold)[1]
    est[delta.idx, "te"]   <- IF.to.est(ifvals$mu_1[delta.idx,] - ifvals$mu_0[delta0.idx,], fold)[1]
    
    se[delta.idx, "se_1"] <- IF.to.est(ifvals$mu_1[delta.idx,] - ifvals$mu_1[delta0.idx,], fold)[2]
    se[delta.idx, "se_0"] <- IF.to.est(ifvals$mu_0[delta.idx,] - ifvals$mu_0[delta0.idx,], fold)[2]
    se[delta.idx, "oe"]   <- IF.to.est(ifvals$mu  [delta.idx,] - ifvals$mu  [delta0.idx,], fold)[2]
    se[delta.idx, "te"]   <- IF.to.est(ifvals$mu_1[delta.idx,] - ifvals$mu_0[delta0.idx,], fold)[2]
    
  }
  
  return(list(est = est, se = se))
  
}
