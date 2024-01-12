setwd("~/research/NPCI/simulation/TPB")

library(dplyr)

### Import help functions
source("Help_function.R")

## Simulation setting
X.trt.vars = c("X1", "X2", "C")
X.out.vars = c("X1", "X2", "C")
rhos = c(0.3, 0.45, 0.6)
rho0 = 0.45

## Generate random sample
data = data.sim(5000, H.nonpar.true, G.nonpar.true)

## Add `g.A` column at data if not already exists
data = data %>%
  group_by(id) %>%
  mutate(g.A = (sum(A) - A) / (n()-1)) %>%
  ungroup()

## Split dat into lists to get N for each cluster
data.list <- split(data, f = data$id)

m <- length(data.list)        # Number of clusters

dummy.mat = matrix(0, nrow = length(rhos), ncol = m)
rownames(dummy.mat) = rhos

estimands  = list(mu   = dummy.mat,
                     mu_1 = dummy.mat,
                     mu_0 = dummy.mat)

for(i in 1:m){
    
  print(paste("     Cluster:", i))
  
  data.i = data.list[[i]]
  
  N = nrow(data.i)
  Y = data.i$Y
  A = data.i$A
  g.A = data.i$g.A
    
  ##### Wrap up internal functions #####
  
  G = function(a) G.nonpar.true(data.i %>% mutate(A = a, g.A = (sum(a)-a)/(length(a)-1)))
  
  pi_ = H.nonpar.true(data.i)
  H = function(a) prod(a*pi_ + (1-a)*(1-pi_))
  
  Q = function(a){
    I(mean(a) >= rho) * H(a) / P.A_bar
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
  
  ##### Replicates over all vectors in A(N) #####
  
  a.rep = expand.grid(replicate(N, 0:1, simplify = F))

  H.rep = apply(X = a.rep, MARGIN = 1, FUN = H)
  
  
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
  
  ### Approximate \Sum_a (w.a + phi.a)^T G.a 
  ### by 1/r * \Sum_{q=1}^{r} (w.a.q + phi.a.q)^T G.a.q / H.a.q
  ### where a.q ~ P(.|X,N) for q = 1,...,r
  
  r = 100
  a_q.rep = t(sapply(1:r,function(q) rbinom(N, 1, pi_)))
  
  ##### Compute IF over rho values #####
  
  for(rho.idx in seq_len(length(rhos))){
    
    rho = rhos[rho.idx]  
    
    print(paste("       rho:", rho))
    
    ### P(A_bar >= rho | X, N)
    P.A_bar = sum(H.rep[which(rowMeans(a.rep) >= rho)])
    
    ##### mu #####
    estimands$mu[rho.idx,i] = 1/N * 1/P.A_bar * sum.part(t=-1)
    
    ##### mu_1 #####
    estimands$mu_1[rho.idx,i] = 1/N * 1/P.A_bar * sum.part(t = 1)
    
    ##### mu_0 #####
    estimands$mu_0[rho.idx,i] = 1/N * 1/P.A_bar * sum.part(t = 0)
    
  }

}

saveRDS(estimands, "estimands.RDS")

