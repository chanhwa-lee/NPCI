#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           Estimation result under CIPS policy             **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********                                                           **********
#**********           Version: 2.0                                    **********
#**********           Jan 12, 2024                                    **********
#*******************************************************************************

############################################
# This file requires "Data/DHS/result.Rdata".
# Causal estimands under the Cluster IPS policy estimation results are visualized.
# Generated figures are saved under "Data/DHS/".
############################################

setwd("~/research/NPCI/application/CIPS")

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tibble)

###------ Read estimator ------###
result.list <- list.files("./Rdata", pattern = "estimate.*rds", full.names = T)
M <- length(result.list)

est.list <- data.frame()
se.list  <- data.frame() 
nuis.list <- list()

for(i in 1:M){
  result = readRDS(result.list[i])
  est.list = rbind(est.list, result$est)
  se.list  = rbind(se.list,  result$se)
  nuis.list = c(nuis.list, result$nuis)
}

print(paste0(M, " estimate Rdata files were loaded"))


###------ Draw plots  ------###

deltas = unique(est.list$delta)

est.med <- se.med <- data.frame(delta = deltas,
                                mu = 0, mu_1 = 0, mu_0 = 0,
                                de = 0, se_1 = 0, se_0 = 0,
                                oe = 0, te = 0)

est.mean <- se.mean <- data.frame(delta = deltas,
                                  mu = 0, mu_1 = 0, mu_0 = 0,
                                  de = 0, se_1 = 0, se_0 = 0,
                                  oe = 0, te = 0)

for(delta.idx in seq_len(length(deltas))){
  
  est.med[delta.idx, ] <- apply(est.list %>% filter(delta == deltas[delta.idx]), 2, median, na.rm = T)
  
  se.med[delta.idx, ] <- sqrt(apply((se.list %>% filter(delta == deltas[delta.idx]))^2, 2, median, na.rm = T))
  
  est.mean[delta.idx, ] <- apply(est.list %>% filter(delta == deltas[delta.idx]), 2, mean, na.rm = T)
  
  se.mean[delta.idx, ] <- sqrt(apply((se.list %>% filter(delta == deltas[delta.idx]))^2, 2, mean, na.rm = T))
  
}


### Drawing arguments ###
limits = c(1/2,2)
breaks = c(1/2,1,2)
mu.limits.y = c(0.675, 0.85)
se1.limits.y = c(-0.06, 0.05)
se0.limits.y = c(-0.06, 0.05)
oe.limits.y  = c(-0.06, 0.05)
de.limits.y  = c(-0.05, 0.125)
te.limits.y  = c(-0.05, 0.125)

### 95% CI plots ###
plot.est <- list()

for(eff in c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")){
  
  data.eff = cbind(est.med %>% select(delta, eff) %>% rename(est.med = eff), 
                   se.med %>% select(eff) %>% rename(se.med = eff)) %>%
    mutate(lower.med = est.med - 1.96 * se.med,
           upper.med = est.med + 1.96 * se.med)
  
  plot.est[[eff]] <- 
    ggplot(data = data.eff) +
    geom_errorbar(aes(x = delta, ymin = lower.med, ymax = upper.med), size = 0.3, width = 0.01) +
    geom_point(aes(x = delta, y = est.med), col = "red", size = 0.5) +
    ggtitle(eff) +
    xlab(expression(delta[0])) +
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
      scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu["CIPS"](delta[0]))),
    plot.est[["mu_1"]] + lims(y = mu.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu["CIPS,1"](delta[0]))),
    plot.est[["mu_0"]] + lims(y = mu.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu["CIPS,0"](delta[0]))),
    plot.est[["de"]] + lims(y = de.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(DE["CIPS"](delta[0]))),
    nrow = 1)
  ,
  plot_grid(
    plot.est[["oe"]] + lims(y = oe.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(OE["CIPS"](delta[0], 1))),
    plot.est[["se_1"]] + lims(y = se1.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(SE["CIPS,1"](delta[0], 1))),
    plot.est[["se_0"]] + lims(y = se0.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(SE["CIPS,0"](delta[0], 1))),
    plot.est[["te"]] + lims(y = te.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(TE["CIPS"](delta[0], 1))),
    nrow = 1)
  ,
  nrow = 2)

ggsave("Fig1.95CIs_CIPS_2x4.pdf", width = 10, height = 6)



###---- D.2. SuperLearner Weights ----###

nuis.data = data.frame()

for(nuis.idx in 1:length(nuis.list)){
  nuis = nuis.list[[nuis.idx]]
  for(split in 1:length(nuis)){
    df = data.frame(nuis[[split]]) %>% rownames_to_column("library")
    df = reshape(df, direction = "long",
                 idvar = "library", ids = "library", 
                 varying = c("pi", "mu"), v.names = "value", 
                 timevar = "type", times = c("pi", "mu"))
    rownames(df) = NULL
    df = df %>% mutate(s = nuis.idx, split = split)
    nuis.data = rbind(nuis.data, df)
  }
}

nuis.data$library <- gsub("^SL\\.|_All$", "", nuis.data$library)
nuis.data$type <- ifelse(nuis.data$type == "mu", "g", nuis.data$type)

nuis.mean.data = nuis.data %>% 
  group_by(s,type,library) %>%
  summarise(value = mean(value))

ggplot(nuis.mean.data, aes(x = library, y = value)) +
  geom_boxplot(fill = "lightblue", color = "steelblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "darkblue", alpha = 0.6) +
  facet_wrap(~ type, scales = "free", labeller = label_parsed) +
  labs(x = "Library", y = "Weight")

ggsave("FigS6.SLWeightCIPS.pdf", width = 8, height = 4)


###---- D.3. Choice of S ----###

S.values = c(5,10,20,30,35,45)

deltas = unique(est.list$delta)

S.max = nrow(est.list)/length(deltas) # 45

est.med <- data.frame()

for(S in S.values){
  
  est.list.temp = est.list[1:(S*length(deltas)),]
  
  est.med.temp <- data.frame(S = S, delta = deltas,
                             mu = 0, mu_1 = 0, mu_0 = 0,
                             de = 0, se_1 = 0, se_0 = 0,
                             oe = 0, te = 0)
  
  for(delta.idx in seq_len(length(deltas))){
    
    est.med.temp[delta.idx, -1] <- apply(est.list.temp %>% filter(delta == deltas[delta.idx]), 2, median, na.rm = T)
    
  }
  
  est.med = rbind(est.med, est.med.temp)
  
}

est.med = est.med %>% select(S, delta, mu, mu_1, mu_0)

est.med.mu = gather(est.med, type, Estimate, mu:mu_0)

est.med.mu$S = as.factor(est.med.mu$S)

mu.limits.y = c(0.75,0.79)
mu1.limits.y = c(0.77,0.81)
mu0.limits.y = c(0.73,0.77)

p.mu = 
  ggplot(data = est.med.mu %>% filter(type == "mu"), 
              aes(x = delta, y = Estimate, group = S, color = S)) +
  geom_line(aes(linetype = S)) +
  xlab(expression(delta[0])) +
  ylab(NULL) + 
  lims(y = mu.limits.y) +
  scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
  ggtitle(expression(mu["CIPS"](delta[0])))

p.mu1 = 
  ggplot(data = est.med.mu %>% filter(type == "mu_1"), 
       aes(x = delta, y = Estimate, group = S, color = S)) +
  geom_line(aes(linetype = S)) +
  xlab(expression(delta[0])) +
  ylab(NULL) + 
  lims(y = mu1.limits.y) +
  scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
  ggtitle(expression(mu["CIPS,1"](delta[0])))

p.mu0 = 
  ggplot(data = est.med.mu %>% filter(type == "mu_0"), 
         aes(x = delta, y = Estimate, group = S, color = S)) +
  geom_line(aes(linetype = S)) +
  xlab(expression(delta[0])) +
  ylab(NULL) + 
  lims(y = mu0.limits.y) +
  scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
  ggtitle(expression(mu["CIPS,0"](delta[0])))


plot_grid(
  p.mu + theme(legend.position="none"),
  p.mu1 + theme(legend.position="none"), 
  p.mu0 + theme(legend.position="none"),
  get_legend(p.mu),
nrow = 1,
rel_widths = c(1,1,1,.3))

ggsave("FigS7.ScompCIPS.pdf", width = 8, height = 4)


