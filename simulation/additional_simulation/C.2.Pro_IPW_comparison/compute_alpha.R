time = proc.time()

###------------------- Load libraries ----------------------###
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))

source("Helpfunc.R")

deltas <- c(0.5,1,2)
###----------- Estimands simulation function ---------------###

m = 10000
N = N.dist(m)
alphas.list <- lapply(N, DGP, type = "alpha.comp", deltas = deltas)
alphas = 
  dplyr::bind_rows(alphas.list) %>%
  group_by(delta) %>%
  summarise(alpha = mean(alpha))

alphas

## Save output
saveRDS(alphas, file = "alphas.rds")


###------------------------------------------------------------###

## Stop timer and report total run time

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))
