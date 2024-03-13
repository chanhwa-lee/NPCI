#*******************************************************************************
#**********           D.4. Estimation result comparison               **********
#**********           with Park et al (2021)                          **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********                                                           **********
#**********           Version: 2.0                                    **********
#**********           Jan 12, 2024                                    **********
#*******************************************************************************

##########################################
# For comparison with the estimation result with Park et al (2021),
# navigate to Chan Park's Github repo https://github.com/qkrcks0218/OptTrt
# and run https://github.com/qkrcks0218/OptTrt/blob/main/Data/9.Summary.R
# est_OMAR.RDS file is generated from 
# by saving `Trt.Prop2` matrix in the script using the following R code
# `saveRDS(Trt.Prop2, file = "est_OMAR.RDS")`
##########################################


setwd("~/research/NPCI/application/CIPS/D.4. Comparison with Park et al (2021)")

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

###------ Read proposed method estimation result ------###
result.list <- list.files("./Rdata", pattern = "estimate.*rds", full.names = T)
M <- length(result.list)

est.list <- data.frame()
se.list  <- data.frame() 
pi.list <- c()

for(i in 1:M){
  result = readRDS(result.list[i])
  est.list = rbind(est.list, result$est)
  se.list  = rbind(se.list,  result$se)
  pi.list = cbind(pi.list, result$pi)
}

print(paste0(M, " estimate Rdata files were loaded"))

est = est.list %>% 
  group_by(delta) %>%
  summarise(mu = median(mu))

se = se.list %>% 
  group_by(delta) %>%
  summarise(mu = median(mu))

pi = apply(pi.list, 1, median)

deltas = est$delta

est_test = est %>% 
  mutate(pi_delta = colMeans(pi %o% deltas / (pi %o% deltas + 1 - pi)),
         se = se$mu,
         mu_l = mu - 1.96*se,
         mu_u = mu + 1.96*se)

saveRDS(est_test, file = "est_CIPS.RDS")


###------ Read from Park et al. (2021) result ------###

est_OMAR = readRDS("est_OMAR.RDS")
est_OMAR = data.frame(est_OMAR[,1:2])
colnames(est_OMAR) = c("mu", "pi")
est_OMAR = est_OMAR %>% mutate(mu = mu/100, type = "OMAR", mu_l = mu, mu_u = mu)
est_OMAR

est_CIPS = est_test %>% 
  filter(pi_delta >= 0.3, pi_delta <= 0.71) %>% 
  mutate(pi = pi_delta,
         type = "CIPS") %>% 
  select(-delta, -pi_delta, -se)

est_CIPS

est_comb = rbind(est_OMAR, est_CIPS)

est_comb

ggplot(data = est_comb %>% filter(pi >= 0.39), aes(x = pi, y = mu, group = type, col = type)) +
  geom_point() +
  geom_line() +
  labs(x = "Average of counterfactual propensity score", y = expression(mu(Q)),
       col = "Policy") +
  scale_color_discrete(labels = c("CIPS", "Park et al")) +
  theme(legend.position = c(0.85, 0.2))

ggsave("FigS9.CIPSvsOMAR.pdf", height = 3, width = 5)
