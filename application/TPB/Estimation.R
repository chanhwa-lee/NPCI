#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           WASH effect estimation                          **********
#**********           under TPB policy                                **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********                                                           **********
#**********           Version: 2.0                                    **********
#**********           Jan 12, 2024                                    **********
#*******************************************************************************

############################################
# This file requires "./Estimator.R" and "../Data/HHData.Rdata".
# Senegal DHS data is analyzed to estimate the causal estimands under TPB policy.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To parallelize, submit jobs in for(s in 1:S){...} separately
# Estimates and SE estimates are saved as "./result/estimate_id?_s?.rds".
############################################

### Load estimator function ###
source("Estimator.R")

### Super Learner libraries ###
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
library(optparse)

SL.library = c("SL.glm", "SL.glmnet", "SL.earth", 
               "SL.gam", "SL.xgboost",
               "SL.ranger", "SL.nnet")
  
### Load preprocessed dataset ###
load("../Data/HHData.Rdata")
data = HH.Data %>% rename(id = cid)

###------------------- Estimator computation ----------------------###

### Parameter setting ###
K = 5
r = 10
S = 1
rhos = 0:60/120
rho0 = 0

time = proc.time()
inter.time = proc.time()

### X.trt.vars: A data.frame with covariates for propensity score model ###
X.trt.vars = colnames(data %>% select(-c(hhid, id, Y, A)))

### X.out.vars: A data.frame with covariates for outcome regression model ###
X.out.vars = colnames(data %>% select(-c(hhid, id, Y, A)))

### Repeat sample splitting and constructing estimators for sample-splitting robust estimator ###
est.list <- list()
se.list <- list()

for(s in 1:S){
  
  print(paste("s:",s))
  
  result <- estimator(data, X.trt.vars, X.out.vars, rhos, rho0, K, r, SL.library)
  
  result = result$estimate

  ### Save output ###
  task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  saveRDS(result, file = paste0("result/estimate_id", task_ID, "_s", s,".rds"))
  
  ### Elapsed time ###
  script.time = proc.time() - inter.time
  inter.time = proc.time()
  print(paste0("Elapsed time: ", script.time[3], " seconds"))
  
}