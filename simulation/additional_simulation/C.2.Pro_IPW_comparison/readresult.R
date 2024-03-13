setwd("~/research/NPCI/simulation/additional_simulation/C.2.Pro_IPW_comparison")
library(dplyr)

###------ Read estimand ------###
file.list <- list.files("../../CIPS/estimand/", pattern = "estimand.*rds", full.names = T)
M <- length(file.list)

deltas <- c(0.5,1,2)
estimands <- data.frame(delta = rep(0, length(deltas)),
                        mu = 0, mu_1 = 0, mu_0 = 0,
                        de = 0, se_1 = 0, se_0 = 0,
                        oe = 0, te = 0)

for(file in file.list){
  estimands <- estimands + readRDS(file)
}

print(paste0(M, " estimand Rdata files were loaded"))

estimands = estimands / M

print(estimands)


###------ Read estimator ------###

### m=500 dataset Proposed & IPW estimators ###

pro.p  = list(est = c(), se = c())
pro.np = list(est = c(), se = c())

ipw.p  = list(est = c(), se = c())
ipw.np = list(est = c(), se = c())

estimate.list = list.files("../../CIPS/data/m500_r100/Rdata/", 
                           pattern = "estimate.*rds", full.names = T)
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


### m=500 dataset Barkley method converged ###

bar.p = list(est = c(), se = c())
bar.cnt = 0

estimate.list = list.files("data/m500_r100/Rdata", 
                           pattern = "estimate.*rds", full.names = T)
M <- length(estimate.list)
for(i in 1:M){
  estimate.i = readRDS(estimate.list[i])
  for(estimate.ii in estimate.i){
    
    if(all(estimate.ii$bar.p$est$mu < 0.1)) next     ## Algorithm not converged
    bar.p$est  = rbind(bar.p$est , estimate.ii$bar.p$est)
    # estimate.ii$bar.p$se = estimate.ii$ipw.np$se %>% 
    #   mutate(mu = 0, mu_1 = 0, mu_0 = 0, de = 0, se_1 = 0, se_0 = 0, oe = 0, te = 0)
    bar.p$se  = rbind(bar.p$se , estimate.ii$bar.p$se)
    bar.cnt = bar.cnt + 1
  }
}

M
bar.cnt

#==============================================================================#


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
      
      # Para result
      apply(est, 2, mean, na.rm = T) - estimand,                                       # Bias
      
      # sqrt(apply((est - rep(1, nrow(est)) %o% estimand)^2, 2, mean, na.rm = T))      # RMSE
      apply(est, 2, sd, na.rm = T),                                                    # SE
      
      colMeans(
        (est - 1.96 * se < rep(1, nrow(est)) %o% estimand) & 
          (rep(1, nrow(est)) %o% estimand < est + 1.96 * se), 
        na.rm = T)                                                                     # 95% Cov
      
    )
    
    result.table = data.frame(result.table)
    
    colnames(result.table) <- c("Truth", "Bias", "ESE", "Cov")
    
    result.table.list = rbind(result.table.list, 
                              data.frame(delta = delta, estimand = rownames(result.table), 
                                         Bias = result.table$Bias, ESE = result.table$ESE, Cov = result.table$Cov)) 
  }
  
  return(result.table.list)
}

data.plot = dplyr::bind_rows(list(pro.p  = result.summarise(pro.p,  estimands),
                                  pro.np = result.summarise(pro.np,  estimands),
                                  ipw.p  = result.summarise(ipw.p,  estimands),
                                  ipw.np = result.summarise(ipw.np,  estimands),
                                  bar.p  = result.summarise(bar.p,   estimands)),
                             .id = "method")

data.plot = data.plot %>% 
  mutate(method = factor(method, 
                         levels = c("pro.np", "ipw.np", "pro.p", "ipw.p", "bar.p"), 
                         ordered = T),
         estimand = factor(estimand, 
                           levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te"),
                           ordered = T))

data.plot.2 = reshape::melt(data.plot, id.vars = c("method", "delta", "estimand"), var = "type")

data.plot.2 = data.plot.2 %>% filter(estimand %in% c("mu", "mu_1", "mu_0"))

### Bar chart (NSS vs P IPW) ###
ggplot(data = data.plot.2 %>% 
         filter(estimand %in% c("mu", "mu_1", "mu_0"), method %in% c("pro.np", "bar.p")) %>%
         mutate(delta = paste0("delta == ", delta)),
       aes(x = estimand, y = value, group = method, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(type ~ delta, scales = "free_y", labeller = label_parsed) +
  scale_x_discrete(labels = c("mu" = expression(mu(delta)),
                              "mu_1" = expression(mu[1](delta)),
                              "mu_0" = expression(mu[0](delta)))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_discrete(name = "Method",
                    labels = c("pro.np" = "NSS",
                               "bar.p"  = "IPW"))

ggsave(filename = "FigS1.CIPS.Pro.Bar.comparison.m500.r100.D1000.pdf", width = 8, height = 5)
