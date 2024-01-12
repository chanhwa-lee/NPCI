time = proc.time()

###------------------- Load libraries ----------------------###
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))


###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-m", "--m"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters for estimand simulation")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters for target estimand computation
r = opt$r              # Number of binary vector sampling
deltas <- c(0.5,1,2)

source("../../Helpfunc.R")

print("[Simulation setting]")
print(paste0("m: ", m))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("r: ", r))

estimands = estimands.sim(m = m, deltas = deltas, r = r)

###-------- Save simulated estimand list as Rdata --------###
## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Save output
saveRDS(estimands, file = paste0("Rdata/estimand_id", task_ID,".rds"))


###------------------------------------------------------------###

## Stop timer and report total run time

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))