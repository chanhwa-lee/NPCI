time = proc.time()

###------------------- Load libraries ----------------------###
library(dplyr)
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(clusteredinterference)
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

suppressPackageStartupMessages(require(optparse))

source("../../Helpfunc.R")

###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-m", "--m"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters per each simulation")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters per simulation
r = opt$r              # Number of binary vector sampling

## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print("[Simulation setting]")
print(paste0("m: ", m))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))

### Simulation ####

D = 1

estimate.save = list()

for(sim.id in 1:D){
  
  print(paste("Sim.id", sim.id))
  
  ## Simulation data generation
  data = data.sim(m)
  print("Data generated")
  
  ## Change format of data to input
  dat <- data %>% 
    select(id, Y, A) %>% 
    group_by(id) %>%
    mutate(g.A = (sum(A) - A) / (n()-1)) %>%
    ungroup()
  
  X.trt <- data %>% select(-c(id, Y, A))
  X.out <- data %>% select(-c(id, Y, A))
  
  estimate.save[[sim.id]] = estimator(dat, X.trt, X.out, deltas, nsplits = 2, r)
  
  print("Estimates and SE estimates computed")
  print("")
}

print(estimate.save)

###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(estimate.save, file = paste0("Rdata/estimate_id", task_ID,".rds"))


## Stop timer and report total run time

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))