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
alphas <- readRDS("../../alphas.rds")

deltas = alphas$delta
alphas = alphas$alpha

## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print("[Simulation setting]")
print(paste0("m: ", m))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("alphas: ", paste(signif(alphas, 4), collapse = ", ")))


### Simulation ####

D = 1

estimate.save = list()

my_grid <- clusteredinterference::makeTargetGrid(alphas = alphas, small_grid = TRUE) 
my_grid = my_grid[1:9,]

for(sim.id in 1:D){
  
  print(paste("Sim.id", sim.id))
  
  ## Simulation data generation
  data = data.sim(m)
  print("Data generated")
  
  ### Proposed & IPW estimators ###
  print("Proposed & IPW estimators")
  
  ## Change format of data to input
  dat <- data %>% 
    select(id, Y, A) %>% 
    group_by(id) %>%
    mutate(g.A = (sum(A) - A) / (n()-1)) %>%
    ungroup()
  
  X.trt <- data %>% select(-c(id, Y, A))
  X.out <- data %>% select(-c(id, Y, A))

  estimate.save[[sim.id]] = estimator(dat, X.trt, X.out, deltas, nsplits = 2, r)
  
  ### Barkley estimators ###
  print("Barkley estimators")
  
  estimate.B <- 
    clusteredinterference::policyFX(data = data, 
                                    formula = Y | A ~ X1 + X2 + (1|id) | id,
                                    target_grid = my_grid, 
                                    k_samps = 1, 
                                    verbose = T)
  
  
  result = cbind(data.frame(delta = rep(deltas,3)),
                 estimate.B$estimates[1:9,c("estimate", "estimand_type")])
  
  estimate.save[[sim.id]]$bar.p$est = 
    tidyr::spread(result, estimand_type, estimate) %>% 
    rename(mu_1 = mu1, mu_0 = mu0) %>%
    select(delta, mu, mu_1, mu_0) %>% 
    mutate(de = 0, se_1 = 0, se_0 = 0, oe = 0, te = 0)
  
  result = cbind(data.frame(delta = rep(deltas,3)),
                 estimate.B$estimates[1:9,c("se", "estimand_type")])

  estimate.save[[sim.id]]$bar.p$se =
    tidyr::spread(result, estimand_type, se) %>%
    rename(mu_1 = mu1, mu_0 = mu0) %>%
    select(delta, mu, mu_1, mu_0) %>%
    mutate(de = 0, se_1 = 0, se_0 = 0, oe = 0, te = 0)
  
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