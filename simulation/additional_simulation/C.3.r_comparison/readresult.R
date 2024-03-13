setwd("~/research/NPCI/simulation/additional_simulation/C.3.r_comparison")

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

m = 100 ## Sample size (number of clusters)

## Help function ##
result.summarise = function(result, estimands){
  
  deltas = estimands$delta
  
  result.table.list = data.frame()
  
  for(delta.idx in seq_len(length(deltas))){
    
    delta = deltas[delta.idx]
    
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
                              result.table %>% 
                                mutate(delta = delta, estimand = rownames(result.table)) %>%
                                select(delta, estimand, Bias, ESE, Cov))
                              
  }
  
  return(result.table.list)
}

data.plot = data.frame()

for(r in c(100, 200, 300, 400, 500)){
  
  est = c()
  se = c()
  
  estimate.list = list.files(paste0("data/m",m,"_r",r,"/Rdata"), 
                             pattern = "estimate.*rds", full.names = T)
  M <- length(estimate.list)
  for(i in 1:M){
    estimate.i = readRDS(estimate.list[i])
    for(estimate.ii in estimate.i){
      est = rbind(est, estimate.ii$est)
      se  = rbind(se , estimate.ii$se)
    }
  }
  
  result = result.summarise(result = list(est = est, se = se), estimands) %>% mutate(r = r)
  
  data.plot = rbind(data.plot, result)

}

data.plot = data.plot %>% 
  mutate(estimand = factor(estimand, 
                           levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te"),
                           ordered = T))

data.plot.2 = reshape::melt(data.plot, id.vars = c("r", "delta", "estimand"), var = "type")

data.plot.2 = data.plot.2 %>% filter(estimand %in% c("mu", "mu_1", "mu_0"))

### Bar chart ###
ggplot(data = data.plot.2 %>% 
         filter(estimand %in% c("mu", "mu_1", "mu_0")) %>%
         mutate(delta = paste0("delta == ", delta)),
       aes(x = estimand, y = value, group = r, fill = r)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(type ~ delta, scales = "free_y", labeller = label_parsed) +
  scale_x_discrete(labels = c("mu" = expression(mu(delta)),
                              "mu_1" = expression(mu[1](delta)),
                              "mu_0" = expression(mu[0](delta)))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(filename = "FigS2.r.comparison.m100.D1000.pdf", width = 8, height = 5)
