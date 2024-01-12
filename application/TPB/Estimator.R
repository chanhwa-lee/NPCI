#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           WASH effect estimation                          **********
#**********           under TPB policy                                **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********                                                           **********
#**********           Version: 2.0                                    **********
#**********           Jan 12, 2024                                    **********
#*******************************************************************************


###----------- Estimator function ---------------###
#' Compute proposed nonparametric estimators from simulated dataset.
#'
#' @param data A data.frame. Columns with cluster `id`, outcome `Y`, treatment `A`, 
#' covariates for treatment model `X.trt.vars`, and covariates for outcome model `X.our.vars`
#' @param X.trt.vars A character vector. Names of treatment model covariates
#' @param X.out.vars A character vector. Names of outcome model covariates
#' @param rhos A numeric vector. rho values of interest for TPB policy
#' @param rho0 A number. A rho value of standard to compute SE, OE, TE
#' @param K An integer. Number of sample splitting groups
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @param SL.library A vector of strings. Names of data-adaptive methods available
#' in `SuperLearner` R package to be used in the ensemble estimation of nuisance
#' functions
#' @return Estimates and se estimates of causal estimands for 6 estimators
#' @examples 
#' data = data.sim(10, H.par.true, G.par.true)
#' X.trt.vars = c("X1", "X2", "C")
#' X.out.vars = c("X1", "X2", "C")
#' rhos = c(0.2, 0.5, 0.8)
#' rho0 = 0.5
#' K = 2
#' r = 10
#' SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet")
#' result = estimator(data, X.trt.vars, X.out.vars, rhos, rho0, K, r, SL.library)
#' result$estimate


estimator <- function(data, X.trt.vars, X.out.vars, rhos, rho0, K, r, SL.library){
  
  ## Record run time
  time = proc.time()
  
  ## Add `g.A` column at data if not already exists
  if(!("g.A" %in% colnames(data))){
    data = data %>%
      group_by(id) %>%
      mutate(g.A = (sum(A) - A) / (n()-1)) %>%
      ungroup()
  }
  
  ## Split dat into lists to get N for each cluster
  data.list <- split(data, f = data$id)
  
  ## setup storage
  n <- length(data.list)        # Number of clusters
  
  ifvals = list(mu   = matrix(0, nrow = length(rhos), ncol = n),
                mu_1 = matrix(0, nrow = length(rhos), ncol = n),
                mu_0 = matrix(0, nrow = length(rhos), ncol = n))
  
  rownames(ifvals$mu) = rhos
  rownames(ifvals$mu_1) = rhos
  rownames(ifvals$mu_0) = rhos
  
  fold <- sample(1:n, replace = F) %% K + 1
  nuis = list()
  
  ###--- Sample Splitting Estimator ---###
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:K) {
    
    print(paste("   Split:", split))
    
    train.idx = which(fold != split)
    eval.idx  = which(fold == split)
    
    data.train = bind_rows(data.list[train.idx])
      
    ### Nonparametric nuisance function estimation (SL)
  
    H.fit.nonpar <- SuperLearner(Y = data.train$A,
                                 X = data.train[,c(X.trt.vars)],
                                 family = binomial(),
                                 SL.library = SL.library)
    
    G.fit.nonpar <- SuperLearner(Y = data.train$Y,
                                   X = data.train[,c(X.out.vars, "A", "g.A")],
                                   family = binomial(),
                                   SL.library = SL.library)
    
    nuis[[split]] = list(pi = H.fit.nonpar$coef, mu = G.fit.nonpar$coef)
    
    ## Evaluate IF values for each cluster in eval fold
    for(i in eval.idx){
      
      print(paste("     Cluster:", i))
      
      data.i = data.list[[i]]
      
      N = nrow(data.i)
      Y = data.i$Y
      A = data.i$A
      g.A = data.i$g.A
      X.out = data.i[,X.out.vars]
      X.trt = data.i[,X.trt.vars]
      
      ##### Wrap up internal functions and      #####
      
      G = function(a){
        newdata = X.out %>% mutate(A = a, g.A = (sum(a) - a) / (length(a) - 1))
        return(predict(G.fit.nonpar, newdata, onlySL = T)$pred[,,drop = T])
      }
      
      pi_ = predict(H.fit.nonpar, X.trt, onlySL = T)$pred[,,drop = T]
      H = function(a) prod(a*pi_ + (1-a)*(1-pi_))
      
      Q = function(a){
        I(mean(a) >= rho) * H.rep[find_row(a)] / P.A_bar
      }
      
      phi.Q = function(a){
        I(mean(a) >= rho) / P.A_bar^2 *
          (all(A == a)*P.A_bar - I(mean(A) >= rho) * H.rep[find_row(a)])
      }
      
      find_row = function(a) sum(a * 2^(seq_along(a)-1)) + 1
      
      w = function(a, t = -1) {
        
        Q.a = Q.rep[find_row(a)]
        
        ### w(a,X,N) for mu(Q)
        if(t == -1){
          return(1/length(a) * rep(Q.a, length(a)))
        }
        
        ### w_t(a,X,N) for mu_t(Q) (t = 0,1)
        Q.a_j = sapply(1:length(a), function(j) {
          a[j] = 1 - a[j]            ## Change value of a_ij to opposite
          Q.rep[find_row(a)]
        })
        return(1 / length(a) * I(a == t) * (Q.a + Q.a_j))
      }
      
      phi = function(a, t = -1){
        
        phi.Q.a = phi.Q.rep[find_row(a)]
        
        ### phi(a,X,N) for mu(Q)
        if(t == -1){
          return(1/length(a) * rep(phi.Q.a, length(a)))
        }
        
        ### phi_t(a,X,N) for mu_t(Q) (t = 0,1)
        phi.Q.a_j = sapply(1:length(a), function(j) {
          a[j] = 1 - a[j]            ## Change value of a_ij to opposite
          phi.Q.rep[find_row(a)]
        })
        return(1 / length(a) * I(a == t) * (phi.Q.a + phi.Q.a_j))
      } 
      
      
      ##### Replicates over all vectors in A(N) #####
      
      a.rep = expand.grid(replicate(N, 0:1, simplify = F))
      
      G.rep = apply(X = a.rep, MARGIN = 1, FUN = G)
      
      H.rep = apply(X = a.rep, MARGIN = 1, FUN = H)
      
      
      ##### Compute IF over rho values #####
      
      for(rho.idx in seq_len(length(rhos))){
        
        rho = rhos[rho.idx]  
        
        ### P(A_bar >= rho | X, N)
        P.A_bar = sum(H.rep[which(rowMeans(a.rep) >= rho)])
        
        ##### Replicates over all vectors in A(N) #####
        
        Q.rep = apply(X = a.rep, MARGIN = 1, FUN = Q)
        
        phi.Q.rep = apply(X = a.rep, MARGIN = 1, FUN = phi.Q)
        
        
        ##### IF of mu #####
        
        ## Compute Outcome Regression part \Sum_a (w.a + phi.a)^T G.a
        w.rep = apply(X = a.rep, MARGIN = 1, FUN = w, t = -1)
        
        phi.rep = apply(X = a.rep, MARGIN = 1, FUN = phi, t = -1)
        
        OR      = sum(colSums((w.rep + phi.rep) * G.rep))
        
        ## Compute Bias Correction part 1/H.hat * w.hat^T (Y - G.hat)
        BC = 1/H(A) * sum(w(A, t = -1) * (Y - G(A)))
        
        ## IF for mu at cluster i
        ifvals$mu[rho.idx, i] <- OR + BC
        
        
        ##### IF of mu_1 #####
        
        ## Compute Outcome Regression part \Sum_a (w_1.a + phi_1.a)^T G.a
        w_1.rep = apply(X = a.rep, MARGIN = 1, FUN = w, t = 1)
        
        phi_1.rep = apply(X = a.rep, MARGIN = 1, FUN = phi, t = 1)
        
        OR_1      = sum(colSums((w_1.rep + phi_1.rep) * G.rep))
        
        ## Compute Bias Correction part 1/H.hat * w_1.hat^T (Y - G.hat)
        BC_1 = 1/H(A) * sum(w(A, t = 1) * (Y - G(A)))
        
        ## IF for mu_1 at cluster i
        ifvals$mu_1[rho.idx, i] <- OR_1 + BC_1
        
        
        ##### IF of mu_0 #####
        
        ## Compute Outcome Regression part \Sum_a (w_0.a + phi_0.a)^T G.a 
        w_0.rep = apply(X = a.rep, MARGIN = 1, FUN = w, t = 0)
        
        phi_0.rep = apply(X = a.rep, MARGIN = 1, FUN = phi, t = 0)
        
        OR_0      = sum(colSums((w_0.rep + phi_0.rep) * G.rep))
        
        ## Compute Bias Correction part 1/H.hat * w_0.hat^T (Y - G.hat)
        BC_0 = 1/H(A) * sum(w(A, t = 0) * (Y - G(A)))
        
        ## IF for mu_0 at cluster i
        ifvals$mu_0[rho.idx, i] <- OR_0 + BC_0
        
      }
      
    }
    
  }
  
  print("")
  print("IF for proposed nonparametric estimator computed.")
  
  ###--- Compute estimates & se from IF ---###
  
  IF.to.est <- function(IF, fold){
    
    est = mean(aggregate(IF, by = list(fold), mean, na.rm = T)[,-1], na.rm = T)
    
    se  = sqrt(1/n * (mean(aggregate(IF, by = list(fold), function(phi) mean(phi^2, na.rm = T))[,-1]) - est^2))
    
    return(c(est, se))
  }
  
  IF.list.to.est <- function(IF.list, fold){
    
    est <- se <-  data.frame(rho = rhos,
                             mu = 0, mu_1 = 0, mu_0 = 0,
                             de = 0, se_1 = 0, se_0 = 0,
                             oe = 0, te = 0)
    
    for(rho.idx in seq_len(length(rhos))){
      
      # rho specific estimates and se estimates
      est[rho.idx, "mu"]   <- IF.to.est(IF.list$mu  [rho.idx,], fold)[1]
      est[rho.idx, "mu_1"] <- IF.to.est(IF.list$mu_1[rho.idx,], fold)[1]
      est[rho.idx, "mu_0"] <- IF.to.est(IF.list$mu_0[rho.idx,], fold)[1]
      est[rho.idx, "de"]   <- IF.to.est(IF.list$mu_1[rho.idx,] - IF.list$mu_0[rho.idx,], fold)[1]
      
      se[rho.idx, "mu"]   <- IF.to.est(IF.list$mu  [rho.idx,], fold)[2]
      se[rho.idx, "mu_1"] <- IF.to.est(IF.list$mu_1[rho.idx,], fold)[2]
      se[rho.idx, "mu_0"] <- IF.to.est(IF.list$mu_0[rho.idx,], fold)[2]
      se[rho.idx, "de"]   <- IF.to.est(IF.list$mu_1[rho.idx,] - IF.list$mu_0[rho.idx,], fold)[2]
      
      # rho versus rho0 comparison estimates and se estimates
      rho0.idx = which(rhos == rho0)
      est[rho.idx, "se_1"] <- IF.to.est(IF.list$mu_1[rho.idx,] - IF.list$mu_1[rho0.idx,], fold)[1]
      est[rho.idx, "se_0"] <- IF.to.est(IF.list$mu_0[rho.idx,] - IF.list$mu_0[rho0.idx,], fold)[1]
      est[rho.idx, "oe"]   <- IF.to.est(IF.list$mu  [rho.idx,] - IF.list$mu  [rho0.idx,], fold)[1]
      est[rho.idx, "te"]   <- IF.to.est(IF.list$mu_1[rho.idx,] - IF.list$mu_0[rho0.idx,], fold)[1]
      
      se[rho.idx, "se_1"] <- IF.to.est(IF.list$mu_1[rho.idx,] - IF.list$mu_1[rho0.idx,], fold)[2]
      se[rho.idx, "se_0"] <- IF.to.est(IF.list$mu_0[rho.idx,] - IF.list$mu_0[rho0.idx,], fold)[2]
      se[rho.idx, "oe"]   <- IF.to.est(IF.list$mu  [rho.idx,] - IF.list$mu  [rho0.idx,], fold)[2]
      se[rho.idx, "te"]   <- IF.to.est(IF.list$mu_1[rho.idx,] - IF.list$mu_0[rho0.idx,], fold)[2]
      
    }
    
    return(list(est = est, se = se))
  }
  
  estimate = IF.list.to.est(ifvals, fold)
  
  result = list(data = data.list,
                rhos = rhos,
                rho0 = rho0,
                fold = fold,
                ifvals = ifvals,
                estimate = estimate,
                nuis = nuis)
  
  script.time = proc.time() - time
  
  print("")
  print("data, fold, IF, estimate, se saved.")
  
  print("")
  print(paste0("Total run time was ", script.time[3], " seconds"))
  
  return(result)
  
}