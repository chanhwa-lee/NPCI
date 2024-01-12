# NPCI
Repository for the paper "**Efficient Nonparametric Estimation of Stochastic Policy Effects with Clustered Interference**"

- Developed nonparametric estimation of causal network effects under interference based on semiparametric efficiency theory and empirical process
- Used ensemble of nonparametric and ML models (spline regression, GAM, boosting, Random Forest, neural net) via SuperLearnear in R
- Assessed the effect of water, sanitation, and hygiene (WASH) facilities on diarrhea incidence among children, allowing for interference within census blocks

## :file_folder: simulation

For simulation, parallel computing at a high-performance computing cluster (Linux, bash) was conducted. 
R files are the main scripts, while bash files (.sh) are for performing parallel computing.

### :file_folder: CIPS

CIPS policy simulation code and results included in the main text

- :page_facing_up: Helpfunc.R: R functions for estimand and estimator computation
- :page_facing_up: estimand.R: Causal estimands approximation
- :page_facing_up: estimator.R: Proposed estimators computation
- :page_facing_up: readresult.R: Read and summarize simulation result
- :file_folder: estimand: Target causal estimands
- :file_folder: data: Estimates of target causal estimands were stored.

### :file_folder: TPB

TPB policy simulation code and results included in the main text

- :page_facing_up: Helpfunc.R: R functions for estimand and estimator computation
- :page_facing_up: estimand.R: Causal estimands approximation
- :page_facing_up: estimator.R: Proposed estimators computation
- :page_facing_up: readresult.R: Read and summarize simulation result
- :page_facing_up: estimands.RDS: Target causal estimands
- :file_folder: data: Estimates of target causal estimands were stored.
  
### :file_folder: additional_simulation

Additional simulation code and results included in the supplementary material

> #### :file_folder: C.1.sizevarydelta
> CIPS with varying delta policy simulation code and results
> - :page_facing_up: Helpfunc.R: R functions for estimand and estimator computation
> - :page_facing_up: estimand.R: Causal estimands approximation
> - :page_facing_up: estimator.R: Proposed estimators computation
> - :page_facing_up: readresult.R: Read and summarize simulation result
> - :file_folder: estimand: Target causal estimands
> - :file_folder: data: Estimates of target causal estimands were stored.

> #### :file_folder: C.2.Pro_IPW_comparison
> Proposed nonparametric estimator versus [Barkley et al. (2020)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-14/issue-3/Causal-inference-from-observational-studies-with-clustered-interference-with-application/10.1214/19-AOAS1314.full)'s IPW estimator comparison
> - :page_facing_up: Helpfunc.R: R functions for estimand and estimator computation
> - :page_facing_up: compute_alpha.R: Compute alpha values in Barkley et al. (2020) corresponding to delta values in CIPS
> - :page_facing_up: alphas.rds: Computed alpha values
> - :page_facing_up: estimator.R: Proposed estimators computation
> - :page_facing_up: readresult.R: Read and summarize simulation result
> - :file_folder: data: Estimates of target causal estimands were stored.
> - NOTE: estimands were not computed here and instead loaded from "~/simulation/CIPS/estimand"

> #### :file_folder: C.3.r_comparison
> Proposed estimator performance over r (subsampling approximation degree) values
> - :page_facing_up: Helpfunc.R: R functions for estimand and estimator computation
> - :page_facing_up: estimator.R: Proposed estimators computation
> - :page_facing_up: readresult.R: Read and summarize simulation result
> - :file_folder: data: Estimates of target causal estimands were stored.
> - NOTE: estimands were not computed here and instead loaded from "~/simulation/CIPS/estimand"

> #### :file_folder: C.4.N_dist_comparison
> Proposed estimator performance over various cluster sizes distributions
> - :file_folder: N3, N3_5, N5, N5_10: Simulation results stored for various cluster sizes distributions. Structure is similar to "~/simulation/CIPS"
> - :page_facing_up: readresult.R: Read and summarize simulation result





## :file_folder: application

Analysis on the effect of water, sanitation, and hygiene (WASH) facilities on diarrhea among children in Senegal accounting for the clustered interference.

### :file_folder: Data

#### Senegal DHS data [(ANSD and ICF, 2020)](https://www.dhsprogram.com/pubs/pdf/FR368/FR368.pdf)
- Sociodemographic, enviromental, and health-related survey on household members 
- Used to assess the effect of WASH facilities on diarrhea incidence among children, allowing for interference within census blocks
- Download the data from [https://dhsprogram.com/data/available-datasets.cfm](https://dhsprogram.com/data/available-datasets.cfm) 
(requires data request submission) and place the datasets in "~/application/Data/" by following procedure:
  - Senegal: Continuous DHS, 2015 -> (download) SNKR7HDT.ZIP -> (uncompress) SNKR7HFL.DTA -> (rename) senegal15.DTA
  - Senegal: Continuous DHS, 2016 -> (download) SNKR7IDT.ZIP -> (uncompress) SNKR7IFL.DTA -> (rename) senegal16.DTA
  - Senegal: Continuous DHS, 2017 -> (download) SNKR7ZDT.ZIP -> (uncompress) SNKR7ZFL.DTA -> (rename) senegal17.DTA
  - Senegal: Continuous DHS, 2018 -> (download) SNKR81DT.ZIP -> (uncompress) SNKR81FL.DTA -> (rename) senegal18.DTA
  - Senegal: Continuous DHS, 2019 -> (download) SNKR8BDT.ZIP -> (uncompress) SNKR8BFL.DTA -> (rename) senegal19.DTA

- :page_facing_up: Preprocessing.R: Preprocessing raw data files to generate HHData.Rdata and generate exploratory figures

### :file_folder: CIPS

CIPS policy application code and results

- :page_facing_up: estimator.R: Proposed estimators computation
- :page_facing_up: Visualization.R: Read and summarize simulation results to generate figures
- :file_folder: Rdata: Estimates of target causal estimands were stored.
- :file_folder: D.4. Comparison with Park et al (2021): Comparison with
  >- :page_facing_up: Preprocessing.R: Preprocessing raw data files to generate HHData.Rdata from "~/application/Data/senegal18.DTA"
  >- :page_facing_up: estimator.R: Proposed estimators computation
  >- :page_facing_up: Visualization.R: Read and summarize simulation results to generate figures
  >- :file_folder: Rdata: Estimates of target causal estimands were stored.

### :file_folder: TPB

TPB policy application code and results

- :page_facing_up: Estimator.R: Proposed estimators computation
- :page_facing_up: Estimation.R: Script for parallelization
- :page_facing_up: Visualization.R: Read and summarize simulation results to generate figures
- :file_folder: result: Estimates of target causal estimands were stored.


- :page_facing_up: Helpfunc.R: R functions for estimand and estimator computation
- :page_facing_up: estimand.R: Causal estimands approximation
- :page_facing_up: estimator.R: Proposed estimators computation
- :page_facing_up: readresult.R: Read and summarize simulation result
- :page_facing_up: estimands.RDS: Target causal estimands
- :file_folder: data: Estimates of target causal estimands were stored.



## :file_folder: Code

### :page_facing_up: Preprocessing.R
- By running this R file, the DHS dataset is cleaned.
- Preprocessed dataset will be saved as "Data/DHS/HHData.Rdata" and used in "Code/Estimation.R".

### :page_facing_up: Estimator.R
- R function for Nonparametric efficient sample splitting estimator under the Cluster Incremental Propensity Score policy.

### :page_facing_up: Estimation.R
- Requires "Code/Estimator.R" and "Data/DHS/HHData.Rdata".
- Senegal DHS data is analyzed to estimate the causal estimands under the Cluster IPS policy.
- The code will take a lot of time, so it is recommended to use parallel computing.
- To parallelize, submit jobs in for(s in 1:S){...} separately.
- Estimates and SE estimates will be saved at "Data/DHS/result.Rdata".

### :page_facing_up: Visualization.R
- Requires "Data/DHS/result.Rdata".
- Causal estimands under the Cluster IPS policy estimation results are visualized.
- Generated figures are saved under "Data/DHS/".


