setwd("~/research/NPCI/simulation/CIPS")
library(dplyr)

###------ Read estimand ------###
file.list <- list.files("estimand", pattern = "estimand.*rds")
M <- length(file.list)

deltas <- c(0.5,1,2)
estimands <- data.frame(delta = rep(0, length(deltas)),
                        mu = 0, mu_1 = 0, mu_0 = 0,
                        de = 0, se_1 = 0, se_0 = 0,
                        oe = 0, te = 0)

for(file in file.list){
  estimands <- estimands + readRDS(paste0("estimand/", file))
}

print(paste0(M, " estimand Rdata files were loaded"))

estimands = estimands / M

print(estimands)

###------ Read estimator ------###

pro.p  = list(est = c(), se = c())
pro.np = list(est = c(), se = c())

ipw.p  = list(est = c(), se = c())
ipw.np = list(est = c(), se = c())

estimate.list = list.files("data/m500_r100/Rdata", pattern = "estimate.*rds", full.names = T)
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


result.summarise = function(result, estimands){
  
  deltas = estimands$delta
  
  result.table.list = data.frame()
  
  for(delta.idx in seq_len(length(deltas))){
    
    delta = deltas[delta.idx]
    # print(paste("delta = ", delta))
    
    est = result$est %>% filter(delta == deltas[delta.idx]) %>% select(-delta)
    
    se = result$se %>% filter(delta == deltas[delta.idx]) %>% select(-delta)
    
    estimand = as.numeric(estimands %>% filter(delta == deltas[delta.idx]) %>% select(-delta))
    
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
                                mutate(delta = delta, estimand = rownames(result.table)) %>%
                                select(estimand, delta, Truth, Bias, RMSE, ASE, ESE, Cov)
                              )
  }
  
  rownames(result.table.list) = NULL
  
  return(result.table.list)
}

result.table = cbind(result.summarise(pro.np,  estimands),
                     result.summarise(pro.p,  estimands) %>% select(Bias:Cov))

result.table$RMSERatio = result.table[,5]/result.table[,10]

result.table = data.frame(lapply(result.table,
                  function(col) if(is.numeric(col)) round(col, 3) else col))

write.csv(result.table, file = "Tab1.CIPS.pro.result.table.m500.r100.D1000.csv")
















