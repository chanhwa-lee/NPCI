time = proc.time()

###------------------- Load libraries ----------------------###
library(dplyr)
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(nnls, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0")
library(gam, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0")
library(SuperLearner, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0")
library(caret)
library(nnet)
library(glmnet)
library(earth)
library(gam)
library(gbm)
library(xgboost)     
library(polspline)
library(ranger)
library(statmod)
library(lme4)

suppressPackageStartupMessages(require(optparse))


###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-n", "--n"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters per each simulation"))
)

opt = parse_args(OptionParser(option_list = option_list))

n = opt$n              # Number of clusters per simulation

### Import help functions
source("../../Help_function.R")

## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Simulation setting
X.trt.vars = c("X1", "X2", "C")
X.out.vars = c("X1", "X2", "C")
rhos = c(0.3, 0.45, 0.6)
rho0 = 0.45
nsplits = 2

print("[Simulation setting]")
print(paste0("n: ", n))
print(paste0("rhos: ", paste(signif(rhos, 4), collapse = ", ")))
print(paste0("K: ", nsplits))


###----------- Simulation Master function ---------------###
#' Compute 6 estimators (proposed, IPW, plug-in) * (Para, Nonpara) from simulated dataset.
#'
#' @param data A data.frame. Columns with cluster `id`, outcome `Y`, treatment `A`, 
#' covariates for treatment model `X.trt.vars`, and covariates for outcome model `X.our.vars`
#' @param X.trt.vars A character vector. Names of treatment model covariates
#' @param X.out.vars A character vector. Names of outcome model covariates
#' @param rhos A numeric vector. rho values of interest for TPB policy
#' @param rho0 A number. A rho value of standard to compute SE, OE, TE
#' @param nsplits A number. Number of sample splits for proposed method
#' @return Estimates and se estimates of causal estimands for 6 estimators
#' @examples 
#' data = data.sim(30, H.par.true, G.par.true)
#' X.trt.vars = c("X1", "X2", "C")
#' X.out.vars = c("X1", "X2", "C")
#' rhos = c(0.2, 0.5, 0.8)
#' rho0 = 0.5
#' nsplits = 2
#' result = simulation.master(data, X.trt.vars, X.out.vars, rhos, rho0, nsplits)
#' result$estimate

simulation.master <- function(data, X.trt.vars, X.out.vars, rhos, rho0, nsplits){
  
  ## Record run time
  time = proc.time()
  
  ## Add `g.A` column at data if not already exists
  data = data %>%
    group_by(id) %>%
    mutate(g.A = (sum(A) - A) / (n()-1)) %>%
    ungroup()
  
  ## Split dat into lists to get N for each cluster
  data.list <- split(data, f = data$id)
  
  ## setup storage
  n <- length(data.list)        # Number of clusters
  N <- sapply(data.list, nrow)  # Number of individuals in each cluster
  
  ifvals = list(mu   = matrix(0, nrow = length(rhos), ncol = n),
                mu_1 = matrix(0, nrow = length(rhos), ncol = n),
                mu_0 = matrix(0, nrow = length(rhos), ncol = n))
  
  rownames(ifvals$mu) = rhos
  rownames(ifvals$mu_1) = rhos
  rownames(ifvals$mu_0) = rhos
  
  ifvals.list = list(pro.p  = ifvals,
                     pro.np = ifvals,
                     ipw.p  = ifvals,
                     ipw.np = ifvals)
  
  fold <- sample(1:n, replace = F) %% nsplits + 1
  
  ###--- Help function to fit nuisance functions ---###
  
  ###--- Help function to compute IF values at each cluster ---###
  eval.IF = function(data.i, G.pred, H.pred){
    
    N = nrow(data.i)
    Y = data.i$Y
    A = data.i$A
    g.A = data.i$g.A
    X.out = data.i[,X.out.vars]
    X.trt = data.i[,X.trt.vars]
    
    ifvals = data.frame(rho = rhos, 
                        mu = rep(0, length(rhos)),
                        mu_1 = rep(0, length(rhos)),
                        mu_0 = rep(0, length(rhos)))
    
    ifvals.list = list(pro = ifvals,
                       ipw = ifvals,
                       plu = ifvals)
    
    ##### Wrap up internal functions and      #####
    
    G = function(a) G.pred(a, X.out)
    
    pi_ = H.pred(X.trt)
    H = function(a) prod(a*pi_ + (1-a)*(1-pi_))
    
    ### rho-specific quantities
    
    Q = function(a){
      I(mean(a) >= rho) * H(a) / P.A_bar
    }
    
    phi.Q = function(a){
      I(mean(a) >= rho) / P.A_bar^2 *
        (all(A == a)*P.A_bar - I(mean(A) >= rho) * H(a))
    }
    
    w = function(a, t = -1) {
      
      Q.a = Q(a)
      
      ### w(a,X,N) for mu(Q)
      if(t == -1){
        return(1/length(a) * rep(Q.a, length(a)))
      }
      
      ### w_t(a,X,N) for mu_t(Q) (t = 0,1)
      Q.a_j = sapply(1:length(a), function(j) {
        a[j] = 1 - a[j]            ## Change value of a_ij to opposite
        Q(a)
      })
      return(1 / length(a) * I(a == t) * (Q.a + Q.a_j))
    }
    
    phi = function(a, t = -1){
      
      phi.Q.a = phi.Q(a)
      
      ### phi(a,X,N) for mu(Q)
      if(t == -1){
        return(1/length(a) * rep(phi.Q.a, length(a)))
      }
      
      ### phi_t(a,X,N) for mu_t(Q) (t = 0,1)
      phi.Q.a_j = sapply(1:length(a), function(j) {
        a[j] = 1 - a[j]            ## Change value of a_ij to opposite
        phi.Q(a)
      })
      return(1 / length(a) * I(a == t) * (phi.Q.a + phi.Q.a_j))
    } 
    
    ##### Functions for approximation #####
    
    kernel.a = function(a,t){
      
      ### Vector [I(1,a_(-j) >= rho) H(1,a_(-j)) + I(0,a_(-j) >= rho) H(0,a_(-j))]_{j=1,...,N} / H(a)
      IH.a.vec = I(mean(a) >= rho) +
        sapply(1:length(a),
               function(j) {
                 a.star = a
                 a.star[j] = 1 - a.star[j]
                 I(mean(a.star) >= rho) * ((1-pi_[j])/pi_[j])^(2*a[j]-1)
               })
      
      return(sum(I(a == t)*IH.a.vec*G(a)))
      
    }
    
    sum.part = function(t = -1){
      if(t == -1){
        return(mean(apply(X = a_q.rep, 
                          MARGIN = 1, 
                          function(a_q) I(mean(a_q) >= rho) * sum(G(a_q)))))
      }
      
      return(mean(apply(X = a_q.rep, 
                        MARGIN = 1, 
                        kernel.a, t)))
    }
    
    diff.part = function(t){
      sum(I(A==t)*G(A)) + 
        sum(
          sapply(1:length(A), function(j){
            A.star = A
            A.star[j] = 1-A[j]
            x = rep(0,N)
            x[j] = I(A.star[j] == t)
            return(sum(x*G(A.star)))
          })
        )
    }
    
    ### Approximate \Sum_a (w.a + phi.a)^T G.a 
    ### by 1/r * \Sum_{q=1}^{r} (w.a.q + phi.a.q)^T G.a.q / H.a.q
    ### where a.q ~ P(.|X,N) for q = 1,...,r
    
    r = 100
    a_q.rep = t(sapply(1:r,function(q) rbinom(N, 1, pi_)))
    
    ##### Compute IF over rho values #####
    
    for(rho.idx in seq_len(length(rhos))){
      
      rho = rhos[rho.idx]  
      
      
      ### P(A_bar >= rho | X, N)
      P.A_bar = mean(sapply(1:10000,function(q) mean(rbinom(N, 1, pi_)) >= rho))
      
      ##### IF of mu #####
      
      ## Compute Outcome Regression part \Sum_a (w.a + phi.a)^T G.a
      
      OR = 
        1/N * 1/P.A_bar * (1 - I(mean(A) >= rho)/P.A_bar) * sum.part(t=-1) + 
        1/N * I(mean(A) >= rho)/P.A_bar * sum(G(A))
      
      ## Compute Bias Correction part 1/H.hat * w.hat^T (Y - G.hat)
      BC = 1/H(A) * sum(w(A, t = -1) * (Y - G(A)))
      
      ## Compute IPW part 1/H.hat * w.hat^T Y
      IPW = 1/H(A) * sum(w(A, t = -1) * Y)
      
      ## IF for mu at cluster i
      ifvals.list$pro$mu[rho.idx] <- OR + BC
      ifvals.list$ipw$mu[rho.idx] <- IPW
      
      
      ##### IF of mu_1 #####
      
      ## Compute Outcome Regression part \Sum_a (w_1.a + phi_1.a)^T G.a
      OR_1 = 
        1/N*1/P.A_bar * (1 - I(mean(A) >= rho)/P.A_bar) * sum.part(t = 1) + 
        1/N*I(mean(A) >= rho)/P.A_bar * diff.part(t = 1)
      
      ## Compute Bias Correction part 1/H.hat * w_1.hat^T (Y - G.hat)
      BC_1 = 1/H(A) * sum(w(A, t = 1) * (Y - G(A)))
      
      ## Compute IPW part 1/H.hat * w_1.hat^T Y
      IPW_1 = 1/H(A) * sum(w(A, t = 1) * Y)
      
      ## IF for mu_1 at cluster i
      ifvals.list$pro$mu_1[rho.idx] <- OR_1 + BC_1
      ifvals.list$ipw$mu_1[rho.idx] <- IPW_1
      
      
      
      ##### IF of mu_0 #####
      
      ## Compute Outcome Regression part \Sum_a (w_0.a + phi_0.a)^T G.a 
      OR_0 = 
        1/N*1/P.A_bar * (1 - I(mean(A) >= rho)/P.A_bar) * sum.part(t = 0) + 
        1/N*I(mean(A) >= rho)/P.A_bar * diff.part(t = 0)
      
      ## Compute Bias Correction part 1/H.hat * w_0.hat^T (Y - G.hat)
      BC_0 = 1/H(A) * sum(w(A, t = 0) * (Y - G(A)))
      
      ## Compute IPW part 1/H.hat * w_0.hat^T Y
      IPW_0 = 1/H(A) * sum(w(A, t = 0) * Y)
      
      ## IF for mu_0 at cluster i
      ifvals.list$pro$mu_0[rho.idx] <- OR_0 + BC_0
      ifvals.list$ipw$mu_0[rho.idx] <- IPW_0
      
    }
    
    return(ifvals.list)
    
  }
  
  
  ###--- Sample Splitting Estimator ---###
  
  print("---- Sample Splitting, IPW, Plug-in Estimator ----")
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:nsplits) {
    
    print(paste("   Split:", split))
    
    train.idx = which(fold != split)
    eval.idx  = which(fold == split)
    
    data.train = bind_rows(data.list[train.idx])
    
    ### 1. Parametric nuisance function estimation
    
      ### 1.(1) H: Main-effect only GLM
      A.formula = as.formula(paste0("A ~" , paste0(X.trt.vars, collapse = "+")))
      H.fit.par <- glm(A.formula, 
                       data = data.train,
                       family = "binomial")
      
      H.pred.par = function(X.trt){
        return(predict(H.fit.par, newdata = X.trt, type = "response"))
      }
      
      ### 1.(2) G: Main-effect only GLM
      Y.formula = as.formula(paste0("Y ~ A + g.A +" , paste0(X.out.vars, collapse = "+")))
      G.fit.par <- glm(Y.formula, 
                       data = data.train, 
                       family = "binomial")
      
      G.pred.par = function(a, X.out){
        newdata = X.out %>% mutate(A = a, g.A = (sum(a) - a) / (length(a) - 1))
        return(predict(G.fit.par, newdata, type = "response"))
      }
      
    ### 2. Nonparametric nuisance function estimation (GMEML, SL)
      
      ### 2.(1) H: GMEML using SL
      H.fit.nonpar <- SuperLearner(Y = data.train$A,
                                   X = data.train[,c(X.trt.vars)],
                                   family = binomial(),
                                   SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
      
      H.pred.nonpar = function(X.trt){
        return(predict(H.fit.nonpar, X.trt, onlySL = T)$pred[,,drop = T])
      }
      
      ### 2.(2) G: SL
      G.fit.nonpar <- SuperLearner(Y = data.train$Y,
                                   X = data.train[,c("A", "g.A", X.out.vars)],
                                   family = binomial(),
                                   SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
      
      G.pred.nonpar = function(a, X.out){
        newdata = X.out %>% mutate(A = a, g.A = (sum(a) - a) / (length(a) - 1))
        return(predict(G.fit.nonpar, newdata, onlySL = T)$pred[,,drop = T])
      }
      
      
    ## Evaluate IF values for each cluster in eval fold
    for(i in eval.idx){
      
      print(paste("     Cluster:", i))
      
      data.i = data.list[[i]]
      
      ### Nonparametric nuisance function estimation IF ###
      
      print(paste("       G: NP, H: NP"))
      
      IF.i.np = eval.IF(data.i, G.pred.nonpar, H.pred.nonpar)
      
      ifvals.list$pro.np$mu[, i]   = IF.i.np$pro$mu
      ifvals.list$pro.np$mu_1[, i] = IF.i.np$pro$mu_1
      ifvals.list$pro.np$mu_0[, i] = IF.i.np$pro$mu_0
      
      ifvals.list$ipw.np$mu[, i]   = IF.i.np$ipw$mu
      ifvals.list$ipw.np$mu_1[, i] = IF.i.np$ipw$mu_1
      ifvals.list$ipw.np$mu_0[, i] = IF.i.np$ipw$mu_0
  
      
      ### Parametric nuisance function estimation IF ###
      
      print(paste("       G: P, H: P"))
      
      IF.i.p = eval.IF(data.i, G.pred.par, H.pred.par)
      
      ifvals.list$pro.p$mu[, i]   = IF.i.p$pro$mu
      ifvals.list$pro.p$mu_1[, i] = IF.i.p$pro$mu_1
      ifvals.list$pro.p$mu_0[, i] = IF.i.p$pro$mu_0
      
      ifvals.list$ipw.p$mu[, i]   = IF.i.p$ipw$mu
      ifvals.list$ipw.p$mu_1[, i] = IF.i.p$ipw$mu_1
      ifvals.list$ipw.p$mu_0[, i] = IF.i.p$ipw$mu_0
      

    }
    
  }
  
  print("")
  print("IF for six estimators (proposed, IPW, plug-in) * (Para, Nonpara) computed.")
  
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
  
  methods = names(ifvals.list)
  estimate =
    lapply(methods,
           function(method) {
             # if (!(method %in% c("pro.p", "pro.np"))) fold = rep(1, n)
             return(IF.list.to.est(ifvals.list[[method]], fold))
           })
  names(estimate) = methods
  
  result = list(data = data.list,
                rhos = rhos,
                rho0 = rho0,
                fold = fold,
                ifvals.list = ifvals.list,
                estimate = estimate)
                # nuis.fit = nuis.fit)
  
  script.time = proc.time() - time
  
  print("")
  print("data, fold, IF, estimate, se, nuisance fit saved.")
  
  print("")
  print(paste0("Total run time was ", script.time[3], " seconds"))
  
  return(result)
  
}





### Simulation ####

m = 1

result.save = list()
estimate.save = list()

for(sim.id in 1:m){
  
  print(paste("Sim.id", sim.id))
  
  ## Simulation data generation
  data = data.sim(n, H.nonpar.true, G.nonpar.true)
  print("Data generated")
  
  ## Run simulation.master function
  result.save[[sim.id]] = simulation.master(data, X.trt.vars, X.out.vars, rhos, rho0, nsplits)
  
  estimate.save[[sim.id]] = result.save[[sim.id]]$estimate
  
  print("Estimates and se estimates computed")
  print("")
}

print(estimate.save)


###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(result.save, file = paste0("Rdata/result_id", task_ID,".rds"))
saveRDS(estimate.save, file = paste0("Rdata/estimate_id", task_ID,".rds"))


###------------------------------------------------------------###

## Stop timer and report total run time

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))