#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           DHS data preprocessing for comparison with      **********
#**********           Park et al (2021). Only using year 2018 data    **********	
#**********                                                           **********	
#**********           Written by:				                              **********
#**********                                                           **********
#**********           Adapted from:                                   **********
#**********           Chan Park(https://github.com/qkrcks0218/OptTrt) **********
#**********                                                           **********
#**********           Version: 2.0                                    **********
#**********           Jan 12, 2024                                    **********
#*******************************************************************************

##########################################
# Download the following data sets from https://dhsprogram.com/data/available-datasets.cfm
# Senegal: Continuous DHS, 2018 -> (download) SNKR81DT.ZIP -> (uncompress) SNKR81FL.DTA -> (rename) senegal18.DTA
# Place the datasets in "~/application/Data/"
##########################################

### Libraries ###
library(readstata13)
library(dplyr)

### Dataset ###
HH.Data.comb = data.frame()

for(year in c(18)){
  print(paste("Year:",year))
  child_filename = paste0("../../Data/senegal", year, ".DTA")
  ChildRaw <- read.dta13(child_filename)
  
  ###### Children data (Age under 60 months) ######
  ChildData <- ChildRaw %>%
    mutate(cid = paste(year, v001, sep = "_"),
           hhid = v002,
           respid = v003,
           C_urban = v025,
           X_resp_age = v012,
           X_HHsize = v136,
           X_HHchildsize = v137,
           X_husb_edu = v701,
           X_husb_job = v705,
           X_resp_edu = v106,
           X_resp_job = v717,
           InterviewDate = v008,
           InterviewMonth = v006,
           X_child_birth = b3,
           X_child_alive = b5,
           X_child_res = b9,
           X_child_age = InterviewDate - X_child_birth,
           X_water = v113,
           X_toilet = v116,
           X_share_t = v160,
           Y_D = h11) %>%
    select(cid, hhid, respid, C_urban, X_resp_age,
           X_HHsize, X_HHchildsize, X_husb_edu, X_husb_job, X_resp_edu, X_resp_job,
           X_child_alive, X_child_res, X_child_age,
           X_water, X_toilet, X_share_t, Y_D)
  
  ChildData = ChildData %>%
    group_by(cid, hhid) %>%
    filter(n_distinct(respid) == 1) %>% ungroup()                               ## Only one mom in the household
  
  ChildData = ChildData %>%                                                     
    na.omit %>%                                                                 ## Complete cases
    filter(X_water != "not de jure",                                            ## remove invalid data
           X_toilet != "not a dejure resident",
           X_share_t != "not a dejure resident",
           X_husb_edu != "don't know",
           X_husb_job != "don't know",
           X_resp_job != "don't know",
           X_child_alive == "yes",
           X_child_res == "respondent")
  
  ###### Outcome ####### 
  ChildData$Y <- 1-as.numeric(ChildData$Y_D=="yes, last two weeks")             ## Indicator whether child had diarrhea in the past two weeks
  
  ###### Treatment ####### 
  ChildData$X_water_ind <- as.numeric(ChildData$X_water %in% 
                                        c("piped into dwelling",
                                          "piped to yard/plot",
                                          "cart with small tank",
                                          "bottled water"))
  
  ChildData$X_share_t <- as.numeric(ChildData$X_share_t=="no")
  
  ChildData$X_toilet <- as.numeric(ChildData$X_toilet %in%
                                     c("flush toilet",
                                       "flush to piped sewer system",
                                       "flush to septic tank",
                                       "flush to pit latrine",
                                       "flush to somewhere else",
                                       "flush, don't know where"))
  
  ChildData$A <- as.numeric( ChildData$X_water_ind 
                             + ChildData$X_toilet*ChildData$X_share_t > 0 )     ## WASH = private water OR flushable toilet
  
  ###### Covariate #######
  ChildData$C_urban <- as.numeric(ChildData$C_urban=="urban")                   ## Indicator of urban area
  
  ChildData$X_resp_edu <- as.numeric(ChildData$X_resp_edu != "no education")    ## Respondent had education
  ChildData$X_husb_edu <- as.numeric(ChildData$X_husb_edu != "no education")    ## Husband had education
  ChildData$X_edu <- I(ChildData$X_resp_edu + ChildData$X_husb_edu > 0)         ## Either one of parents had education
  
  ChildData$X_resp_job <- as.numeric(ChildData$X_resp_job != "not working")     ## Respondent has job
  ChildData$X_husb_job <- as.numeric(ChildData$X_husb_job != "did not work")    ## Husband has job
  ChildData$X_job <- I(ChildData$X_resp_job + ChildData$X_husb_job > 0)         ## Either one of parents has job                                  
  
  Data = ChildData
  
  ###### Household-level Data #######
  HH.Data <- Data %>%                                                           ## Average child-level covariates within a household
    select(cid, hhid, Y, A, 
           X_HHsize, X_HHchildsize, 
           X_job, X_edu, X_resp_age, X_child_age, C_urban) %>%
    group_by(cid, hhid) %>% 
    summarise(Y = floor(mean(Y)),
              A = mean(A),
              X_HHsize = mean(X_HHsize),
              X_HHchildsize = mean(X_HHchildsize),
              X_HHchildsizereal = n(),
              X_job = mean(X_job),
              X_edu = mean(X_edu),
              X_resp_age = mean(X_resp_age),
              X_child_age = mean(X_child_age),
              C_urban = mean(C_urban)) %>%
    filter(X_HHchildsize == X_HHchildsizereal) %>%
    group_by(cid) %>%
    mutate(C_clsize = n()) %>%
    ungroup() %>%
    filter(C_clsize > 1) %>%
    select(-X_HHchildsizereal) %>%
    na.omit
  
  ###### Summary statistics #######
  print(paste("Y prop =", round(mean(HH.Data$Y),3)))
  print(paste("A prop =", round(mean(HH.Data$A),3)))
  print(paste("Corr =", round(cor(HH.Data$Y, HH.Data$A),3)))
  
  print(summary(glm(A ~ . -cid -hhid -Y, 
                    data = HH.Data, family = "binomial")))                        
  
  print(summary(glm(Y ~ . -cid -hhid, 
                    data = HH.Data, family = "binomial")))
  
  ###### Combine Household-level Datasets #######
  HH.Data.comb = rbind(HH.Data.comb, HH.Data)
}

### Year-combined analysis ###
HH.Data = HH.Data.comb
print(summary(glm(A ~ . -cid -hhid -Y, 
                  data = HH.Data, family = "binomial")))                        

print(summary(glm(Y ~ . -cid -hhid, 
                  data = HH.Data, family = "binomial")))

print(summary(glm(Y ~ . -cid -hhid, 
                  data = HH.Data %>% group_by(cid) %>%
                    mutate(g.A = ifelse(n() == 1, 0, (sum(A) - A) / (n()-1))) %>%
                    ungroup(), family = "binomial")))


length(unique(HH.Data$cid))                                                     ## 1074 clusters
table(table(HH.Data$cid))                                                       ## Distribution of cluster size 
hist(table(HH.Data$cid), breaks = 0:12)                                         ## Distribution cluster size

### Save preprocessed dataset ###
save(HH.Data, file = "HHData.Rdata")
