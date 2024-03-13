#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           Estimation result under TPB policy              **********
#**********           visualization                                   **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********                                                           **********
#**********           Version: 2.0                                    **********
#**********           Jan 12, 2024                                    **********
#*******************************************************************************


library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tibble)

###------ Read estimator ------###
result.list <- list.files("./result", pattern = "estimate.*rds", full.names = T)
M <- length(result.list)

est.list <- data.frame()
se.list  <- data.frame() 

for(i in 1:M){
  result = readRDS(result.list[i])
  est.list = rbind(est.list, result$est)
  se.list  = rbind(se.list,  result$se)
}

print(paste0(M, " estimate Rdata files were loaded"))

est.list.comb = bind_rows(est.list)
se.list.comb  = bind_rows(se.list)

rhos = 0:60/120
rho0 = 0

est.med <- se.med <- data.frame(rho = rhos,
                                mu = 0, mu_1 = 0, mu_0 = 0,
                                de = 0, se_1 = 0, se_0 = 0,
                                oe = 0, te = 0)

median((est.list.comb %>% filter(rho == 0))$mu_0)
median((est.list.comb %>% filter(rho == 1/120))$mu_0)

for(rho.idx in seq_len(length(rhos))){

  est.med[rho.idx, ] <- apply(est.list.comb %>% filter(rho == rhos[rho.idx]), 2, median)
  se.med[rho.idx, ]  <- apply(se.list.comb  %>% filter(rho == rhos[rho.idx]), 2, median)

}

result = list(est.list = est.list, se.list = se.list,
              est.med = est.med, se.med = se.med)


### Estimation result ###
est.med = result$est.med
se.med = result$se.med

### Drawing arguments ###
limits = c(0,0.5)
breaks = 0:10/10
mu.limits.y = c(0.70, 0.85)
oe.limits.y = c(-0.025, 0.03)
se1.limits.y = c(-0.025, 0.03)
se0.limits.y = c(-0.025, 0.03)
de.limits.y = c(0, 0.10)
te.limits.y = c(0, 0.10)

### 95% CI plots ###
plot.est <- list()

for(eff in c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")){
  
  data.eff = cbind(est.med %>% select(rho, eff) %>% rename(est.med = eff), 
                   se.med %>% select(eff) %>% rename(se.med = eff)) %>%
    mutate(lower.med = est.med - 1.96 * se.med,
           upper.med = est.med + 1.96 * se.med)
  
  plot.est[[eff]] <- 
    ggplot(data = data.eff) +
    geom_errorbar(aes(x = rho, ymin = lower.med, ymax = upper.med), size = 0.3, width = 0.005) +
    geom_point(aes(x = rho, y = est.med), col = "red", size = 0.5) +
    ggtitle(eff) +
    xlab(expression(rho)) +
    ylab(NULL)
  
  if(eff %in% c("de", "se_1", "se_0", "oe", "te")){
    plot.est[[eff]] <- plot.est[[eff]] + geom_hline(yintercept = 0, color = "blue", size = 0.5, linetype = "dashed")
  }
  
}

### Arranged plots (all 8 plots: 2x4) ###
plot_grid(
  plot_grid(
    plot.est[["mu"]] + 
      lims(y = mu.limits.y) + 
      scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(mu["TPB"](rho))),
    plot.est[["mu_1"]] + lims(y = mu.limits.y) + scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(mu["TPB,1"](rho))),
    plot.est[["mu_0"]] + lims(y = mu.limits.y) + scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(mu["TPB,0"](rho))),
    plot.est[["de"]] + lims(y = de.limits.y) + scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(DE["TPB"](rho))),
    nrow = 1)
  ,
  plot_grid(
    plot.est[["oe"]] + lims(y = oe.limits.y) + scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(OE["TPB"](rho, 0))),
    plot.est[["se_1"]] + lims(y = se1.limits.y) + scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(SE["TPB,1"](rho, 0))),
    plot.est[["se_0"]] + lims(y = se0.limits.y) + scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(SE["TPB,0"](rho, 0))),
    plot.est[["te"]] + lims(y = te.limits.y) + scale_x_continuous( limits = limits, breaks = breaks) + 
      ggtitle(expression(TE["TPB"](rho, 0))),
    nrow = 1),
  nrow = 2)

ggsave("Fig2.95CIs_TPB_2x4.pdf", width = 10, height = 6)





###---- D.3. Choice of S ----###

S.values = c(5,10,20,30,35,45)

rhos = unique(est.list$rho)

S.max = nrow(est.list)/length(rhos) # 50

est.med <- data.frame()

for(S in S.values){
  
  est.list.temp = est.list[1:(S*length(rhos)),]
  
  est.med.temp <- data.frame(S = S, rho = rhos,
                             mu = 0, mu_1 = 0, mu_0 = 0,
                             de = 0, se_1 = 0, se_0 = 0,
                             oe = 0, te = 0)
  
  for(rho.idx in seq_len(length(rhos))){
    
    est.med.temp[rho.idx, -1] <- apply(est.list.temp %>% filter(rho == rhos[rho.idx]), 2, median, na.rm = T)
    
  }
  
  est.med = rbind(est.med, est.med.temp)
  
}

est.med = est.med %>% select(S, rho, mu, mu_1, mu_0)

est.med.mu = gather(est.med, type, Estimate, mu:mu_0)

est.med.mu$S = as.factor(est.med.mu$S)

mu.limits.y = c(0.775,0.795)
mu1.limits.y = c(0.795,0.815)
mu0.limits.y = c(0.745,0.765)

p.mu = 
  ggplot(data = est.med.mu %>% filter(type == "mu"), 
         aes(x = rho, y = Estimate, group = S, color = S)) +
  geom_line(aes(linetype = S)) +
  xlab(expression(rho)) +
  ylab(NULL) + 
  lims(y = mu.limits.y) +
  scale_x_continuous(limits = limits, breaks = breaks) + 
  ggtitle(expression(mu["TPB"](rho)))

p.mu1 = 
  ggplot(data = est.med.mu %>% filter(type == "mu_1"), 
         aes(x = rho, y = Estimate, group = S, color = S)) +
  geom_line(aes(linetype = S)) +
  xlab(expression(rho)) +
  ylab(NULL) + 
  lims(y = mu1.limits.y) +
  scale_x_continuous(limits = limits, breaks = breaks) + 
  ggtitle(expression(mu["TPB,1"](rho)))

p.mu0 = 
  ggplot(data = est.med.mu %>% filter(type == "mu_0"), 
         aes(x = rho, y = Estimate, group = S, color = S)) +
  geom_line(aes(linetype = S)) +
  xlab(expression(rho)) +
  ylab(NULL) + 
  lims(y = mu0.limits.y) +
  scale_x_continuous(limits = limits, breaks = breaks) + 
  ggtitle(expression(mu["TPB,0"](rho)))


plot_grid(
  p.mu + theme(legend.position="none"),
  p.mu1 + theme(legend.position="none"), 
  p.mu0 + theme(legend.position="none"),
  get_legend(p.mu),
  nrow = 1,
  rel_widths = c(1,1,1,.3))

ggsave("FigS8.ScompTPB.pdf", width = 8, height = 4)
