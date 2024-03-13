setwd("~/research/NPCI/simulation/TPB")
library(dplyr)

rhos = c(0.3, 0.45, 0.6)
rho0 = 0.45

estimands.raw = readRDS("estimands.RDS")

estimands <- data.frame(rho = rhos,
                        mu = 0, mu_1 = 0, mu_0 = 0)

estimands$mu   = rowMeans(estimands.raw$mu)
estimands$mu_1 = rowMeans(estimands.raw$mu_1)
estimands$mu_0 = rowMeans(estimands.raw$mu_0)

estimands = estimands %>% 
  mutate(de = mu_1 - mu_0,
         se_1 = mu_1 - estimands[estimands$rho == rho0,"mu_1"],
         se_0 = mu_0 - estimands[estimands$rho == rho0,"mu_0"],
         oe = mu - estimands[estimands$rho == rho0,"mu"],
         te = mu_1 - estimands[estimands$rho == rho0,"mu_0"])

estimands

###------ Read estimator ------###

pro.p  = list(est = c(), se = c())
pro.np = list(est = c(), se = c())

ipw.p  = list(est = c(), se = c())
ipw.np = list(est = c(), se = c())

estimate.list = list.files("data/n500/Rdata", pattern = "estimate.*rds", full.names = T)
M <- length(estimate.list)
for(i in 1:M){
  estimate.i = readRDS(estimate.list[i])
  for(estimate.ii in estimate.i){
    pro.p$est  = rbind(pro.p$est , estimate.ii$pro.p$est)
    pro.p$se  = rbind(pro.p$se , estimate.ii$pro.p$se)
    
    pro.np$est = rbind(pro.np$est, estimate.ii$pro.np$est)
    pro.np$se = rbind(pro.np$se, estimate.ii$pro.np$se)
    
    ipw.p$est  = rbind(ipw.p$est , estimate.ii$ipw.p$est)
    ipw.p$se  = rbind(ipw.p$se , estimate.ii$ipw.p$se)
    
    ipw.np$est = rbind(ipw.np$est, estimate.ii$ipw.np$est)
    ipw.np$se = rbind(ipw.np$se, estimate.ii$ipw.np$se)
    
  }
}


###------ Result table ------###

result.summarise = function(result, estimands){
  
  rhos = estimands$rho
  
  result.table.list = data.frame()
  
  for(rho.idx in seq_len(length(rhos))){
    
    rho = rhos[rho.idx]
    # print(paste("rho = ", rho))
    
    est = result$est %>% filter(rho == rhos[rho.idx]) %>% select(-rho)
    
    se = result$se %>% filter(rho == rhos[rho.idx]) %>% select(-rho)
    
    estimand = as.numeric(estimands %>% filter(rho == rhos[rho.idx]) %>% select(-rho))
    
    result.table <- cbind(
      
      # True Estimand
      estimand,                                                                    
      
      apply(est, 2, mean) - estimand,                                          # Bias
      
      sqrt(apply((est - rep(1, nrow(est)) %o% estimand)^2, 2, mean)),      # RMSE
      
      apply(se, 2, mean)   ,                                                   # Avg. SE
      
      apply(est, 2, sd)    ,                                                   # Emp. SE
      
      colMeans(
        (est - 1.96 * se < rep(1, nrow(est)) %o% estimand) & 
          (rep(1, nrow(est)) %o% estimand < est + 1.96 * se), 
        na.rm = T)                                                                  # 95% coverage
      
    )
    
    result.table = data.frame(result.table)
    
    colnames(result.table) <- c("Truth", "Bias", "RMSE", "ASE", "ESE", "Cov")
    
    result.table.list = rbind(result.table.list, 
                              result.table %>%
                                mutate(rho = rho, estimand = rownames(result.table)) %>%
                                select(estimand, rho, Truth, Bias, RMSE, ASE, ESE, Cov)
    )
  }
  
  rownames(result.table.list) = NULL
  
  return(result.table.list)
}


##-- Proposed NP vs P --##

result.table = cbind(result.summarise(pro.np,  estimands),
                     result.summarise(pro.p,  estimands) %>% select(Bias:Cov))

result.table$RMSERatio = result.table[,5]/result.table[,10]

result.table = data.frame(lapply(result.table,
                                 function(col) if(is.numeric(col)) round(col, 3) else col))

result.table

write.csv(result.table, file = "Tab2.TPB.pro.result.table.m500.r100.D1000.csv")

