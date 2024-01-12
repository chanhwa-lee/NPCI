### Parameter setting ###
nsplits = 5
r = 5000
s = 1
delta_l = 1/2
delta_u = 2
num.delta = 30

library(dplyr)
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(nnls, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0")
library(SuperLearner, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0")
library(ggplot2, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0")
library(SuperLearner)
library(ggplot2)
library(caret)
library(nnet)
library(glmnet)
library(earth)
library(gam)
library(gbm)
library(xgboost)     
library(polspline)
library(ranger)

SL.library = c("SL.glm", "SL.glmnet", "SL.earth", 
               "SL.gam", "SL.xgboost",
               "SL.ranger", "SL.nnet")

### Dataset ###
load("../Data/HHData.Rdata")

###------------------ Estimator function -------------------###

estimator <- function(dat, X.trt, X.out, deltas, nsplits, r, SL.library){
  
  ## If deltas does not include 1, then add it to the list
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  
  ## Split dat into lists to get N for each cluster
  dat.list <- split(dat, f = dat$id)
  
  ## setup storage
  n <- length(dat.list)        # Number of clusters
  N <- sapply(dat.list, nrow)  # Number of individuals in each cluster
  N.limit <- cumsum(c(0,N))
  
  est <- se <-  data.frame(delta = deltas,
                           mu = 0, mu_1 = 0, mu_0 = 0,
                           de = 0, se_1 = 0, se_0 = 0,
                           oe = 0, te = 0)
  
  pi <- numeric(sum(N))
  mu <- numeric(sum(N))
  
  ifvals.mu   <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu per cluster
  ifvals.mu_1 <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu_1 per cluster
  ifvals.mu_0 <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu_0 per cluster
  
  fold <- sample(1:n, replace = F) %% nsplits + 1
  foldlong <- rep(fold, times = N)
  
  nuis = list()
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:nsplits) {
    
    print(paste("   Split:", split))
    
    train.idx = which(foldlong != split)
    eval.idx  = which(foldlong == split)
    
    ### 2. SuperLearner (RF & GLM)
    pi.fit = SuperLearner(Y = dat$A[train.idx], 
                          X = X.trt[train.idx,],
                          family = binomial(), 
                          SL.library = SL.library)
    
    pi[eval.idx] <- predict(pi.fit, X.trt[eval.idx,], onlySL = TRUE)$pred
    
    mu.fit <- SuperLearner(Y = dat$Y[train.idx], 
                           X = cbind(A = dat$A, g.A = dat$g.A, X.out)[train.idx,],
                           family = binomial(), 
                           SL.library = SL.library)
    
    mu[eval.idx] <- predict(mu.fit, cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], onlySL = TRUE)$pred
    
    nuis[[split]] = list(pi = pi.fit$coef, mu = mu.fit$coef)
    
    ## evaluate IF values for each cluster in eval fold
    for(i in which(fold == split)){
      
      # print(paste("   Cluster:", i))
      
      dat.sub.idx <- (N.limit[i]+1):N.limit[i+1]
      dat.sub <- dat[dat.sub.idx,]
      
      pi.sub  <- pi[dat.sub.idx]
      mu.sub  <- mu[dat.sub.idx]
      
      if(2^N[i] <= r){
        
        ### Summing all vectors in A(N)
        a = as.vector(t(expand.grid(replicate(N[i], 0:1, simplify = F))))
        rr = 2^N[i]                                                             ## rr: placeholder for r
        
      }else{
        
        ### Subsampling approximation: Random binary vector from A(N_i) w. uniform sampling
        a <- sample(c(0,1), r*N[i], replace = T)                                # rN[i] vector
        rr = r                                                                  ## rr: placeholder for r
        
      }
      
      # Help function
      agg.fun <- function(a, fun){
        aggregate(a, list(rep(1:rr, each = N[i])), fun)[,-1]
      }
      
      g.a = (rep(agg.fun(a, sum), each = N[i]) - a) / (N[i]-1)
      
      ## Outcome regression
      mu.a <- predict(mu.fit, 
                      cbind(data.frame(A = a, g.A = g.a), bind_rows(replicate(rr, X.out[dat.sub.idx,], simplify = FALSE))),
                      onlySL = TRUE)$pred
      
      
      pi.delta.sub <- pi.sub %o% deltas / (pi.sub %o% deltas + 1 - pi.sub)       # N[i] x length(deltas) matrix
      
      pi.delta.sub.rep <- rep(1, rr) %x% pi.delta.sub                             # rN[i] x length(deltas) matrix (replicate pi.delta.sub r times)
      
      cluster.prob <- agg.fun(a*pi.delta.sub.rep + (1-a)*(1-pi.delta.sub.rep), prod)   # r x length(deltas) matrix
      
      add.factor <-
        (2*a-1) %o% deltas * rep(dat.sub$A - pi.sub, rr) / 
        (rep(1, rr) %x% (pi.sub %o% deltas + 1 - pi.sub)) /
        (a*rep(pi.sub, rr) %o% deltas + (1-a)*(1-rep(pi.sub, rr)))                 # rN[i] x length(deltas) matrix
      
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
      
      ## Bias correction
      cluster.prob.ratio <- deltas^(sum(dat.sub$A)) / 
        apply(pi.sub %o% deltas + 1 - pi.sub, 2, prod)                           # length(deltas) vector
      
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
      ifvals.mu[,i]   <- out.reg.mu + bias.cor.mu
      
      ## IF for mu_1 at cluster i
      ifvals.mu_1[,i] <- out.reg.mu_1 + bias.cor.mu_1
      
      ## IF for mu_0 at cluster i
      ifvals.mu_0[,i] <- out.reg.mu_0 + bias.cor.mu_0
      
    }
    
  }
  
  ## helper function
  IF.to.est <- function(IF){
    est = mean(aggregate(IF, by = list(fold), mean)[,-1])
    se  = sqrt(1/n * (mean(aggregate(IF, by = list(fold), function(phi) mean(phi^2))[,-1]) - est^2))
    return(c(est, se))
  }
  
  ## Standard (delta = 1)
  delta1.idx = which(deltas == 1)
  
  for(delta.idx in seq_len(length(deltas))){
    
    
    ### Nonparametric
    # delta specific estimates and se estimates
    est[delta.idx, "mu"]   <- IF.to.est(ifvals.mu  [delta.idx,])[1]
    est[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1[delta.idx,])[1]
    est[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0[delta.idx,])[1]
    est[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta.idx,])[1]
    
    se[delta.idx, "mu"]   <- IF.to.est(ifvals.mu  [delta.idx,])[2]
    se[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1[delta.idx,])[2]
    se[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0[delta.idx,])[2]
    se[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta.idx,])[2]
    
    # delta versus delta==1 comparison estimates and se estimates
    est[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_1[delta1.idx,])[1]
    est[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0[delta.idx,] - ifvals.mu_0[delta1.idx,])[1]
    est[delta.idx, "oe"]   <- IF.to.est(ifvals.mu  [delta.idx,] - ifvals.mu  [delta1.idx,])[1]
    est[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta1.idx,])[1]
    
    se[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_1[delta1.idx,])[2]
    se[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0[delta.idx,] - ifvals.mu_0[delta1.idx,])[2]
    se[delta.idx, "oe"]   <- IF.to.est(ifvals.mu  [delta.idx,] - ifvals.mu  [delta1.idx,])[2]
    se[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta1.idx,])[2]
    
  }
  
  #  return(list(dat = dat, X.trt = X.trt, X.out = X.out, fold = fold, estimates = estimates, var = var))
  return(list(est = est, se = se, nuis = nuis))
  
}

###------------------- Estimator computation ----------------------###

time = proc.time()
inter.time = proc.time()

dat <- HH.Data %>% 
  group_by(cid) %>%
  mutate(g.A = ifelse(n() == 1, 0, (sum(A) - A) / (n()-1))) %>%
  ungroup() %>%
  select(cid, Y, A, g.A) %>% 
  rename(id = cid)

X.trt <- HH.Data %>%
  select(-c(hhid, cid, Y, A))

X.out <- HH.Data %>%
  select(-c(hhid, cid, Y, A))

deltas <- c(exp(seq(from = log(delta_l), to = 0, length.out = num.delta))[-num.delta],
            1,
            exp(seq(from = 0, to = log(delta_u), length.out = num.delta))[-1])

## Repeat computing estimator for robustness from sample splitting
est.list <- list()
se.list <- list()
nuis.list <- list()

for(k in 1:s){
  
  print(paste("s:",k))
  
  result <- estimator(dat, X.trt, X.out, deltas, nsplits, r, SL.library)
  
  est.list[[k]] <- result$est
  se.list[[k]] <- result$se
  nuis.list[[k]] <- result$nuis
  
  script.time = proc.time() - inter.time
  inter.time = proc.time()
  print(paste0("Elapsed time: ", script.time[3], " seconds"))
  
}

result = list(est = bind_rows(est.list), 
              se  = bind_rows(se.list),
              nuis = nuis.list)

### Save output ###
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
saveRDS(result, file = paste0("Rdata/estimate_id", task_ID,".rds"))

### Elapsed time ###
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))