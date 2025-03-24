# Script to run the brms model for the SpaNov paper (novelyScore)
# Date:  11/10/2023
# Author: Joern Alexander Quent
# /* 
# ----------------------------- Need to be changed ---------------------------
# */

# Use correct location based on computerS
if(Sys.info()[4] == "DESKTOP-335I26I"){
  # Work laptop
  setwd("D:/researchProjects/OLM_project/analysis")
} else if(Sys.info()[4] == 'alex-Zenbook-UX3404VA-UX3404VA') {
   setwd("/media/alex/shared/researchProjects/OLM_project/analysis")
}  else {
  #setwd("~/GitHub/OLM_project/analysis") # for work desktop
  setwd("~/Work/researchProjects/OLM_project/analysis") # for my laptop
  #setwd("C:/Users/Gaming machine/Desktop/Alex Work/OLM_project/analysis")  # home desktop
}


# /* 
# ----------------------------- Notes on this script ---------------------------
# */

# /* 
# ----------------------------- BRMS parameters ---------------------------
# */
iter = 10000
chains = 10
cores = 10

# /* 
# ----------------------------- Libraries & functions ---------------------------
# */
library(brms)
library(ggplot2)
library(plyr)
library(assortedRFunctions)
library(stringr)
library(data.table)

# /* 
# ----------------------------- Load data and prepare ---------------------------
# */
# Load 
load("ignore_eventTable3/images/SectorSize_10_SpaNov_contPM_smo6.RData")

# Load other data
path2data <- "data/ignore_fMRI_version1/"
load(paste0(path2data, "combined_data/OLM_7T_all_data.RData"))
lookupTable  <- read.csv(paste0(path2data, "lookUpTable.csv"))

# Load the subjects that are included in this analysis
subjectFile  <- readLines(paste0(path2data, "OHBM_subject2analyse.txt"))
subjIDs_R    <- str_split(subjectFile, pattern = ",")[[1]] 

# Select the subject IDs to include in this analysis
subjIDs <- lookupTable$anonKey[lookupTable$Rnum %in% subjIDs_R]

# Subset to data that is being included in the analysis
OLM_7T_trial_results <- OLM_7T_trial_results[OLM_7T_trial_results$subject %in% subjIDs, ]

# Convert list to data frame
data <- rbindlist(subj_list, idcol = "subject")

# Add subjects' R number
data$ppid <- subjIDs_R[data$subject]

# Function to match runStartTime
find_runStartTime <- function(ppid, run){
  # Get corresponding to find the run in the trial data
  anonKey <- lookupTable$anonKey[lookupTable$Rnum == ppid[1]]
  
  # Use the anonKey & run to get runStartTime
  runStartTime <- OLM_7T_trial_results$runStartTime[OLM_7T_trial_results$subject == anonKey & 
                                                    OLM_7T_trial_results$block_num == run[1]]
  
  return(runStartTime[1])
}

data <- ddply(data, c("subject", "ppid", "run"), mutate, 
              runStartTime = find_runStartTime(ppid, run))

data$subject <- as.character(data$subject)

# Make time relative to the start of the run. The real onset times will be 
# slightly different but this will not matter for this. 
data$onset_rel <- data$onset - data$runStartTime

# Add run type
data$runType <- "encoding"
data$runType[data$run == 2 | data$run == 4] <- "retrieval"

# Subset to only encoding
data_sub <- data[data$runType == "encoding", ]

# Change run number to match the description in paper. Run 3 is Encoding Run 2
data_sub$run[data_sub$run == 3] <- 2

# Convert run to factor
data_sub$f_run <- as.factor(data_sub$run)

# /* 
# ----------------------------- OLM_7T novelty score ~ time since onset for Run 1 & Run 2 ---------------------------
# */
# Set brms priors
prior_noveltyScore  <- c(prior(normal(0, 1), class = "Intercept"),
                         prior(normal(0, 1), class = "b"))

# Run 1 
## Remove NA and get correct run
data_sub3 <- na.omit(data_sub[data_sub$run == 1, ])

## Scale x and y
data_sub3$s_onset_rel    <- scale_m0_sd0.5(data_sub3$onset_rel)
data_sub3$s_noveltyScore <- scale(data_sub3$noveltyScore)

## Fit brms model
m_noveltyScore_run1 <- brm(bf(s_noveltyScore  ~ s_onset_rel  + (s_onset_rel | subject)), 
          data = data_sub3,
          family = student(), 
          prior = prior_noveltyScore, 
          iter = iter,
          chains = chains,
          cores = cores,
          seed = 1125,
          control = list(adapt_delta = .95),
          silent = 0,
          refresh = 1,
          backend = "cmdstanr")

summary(m_noveltyScore_run1)

# Run 2 
## Remove NA and get correct run
data_sub3 <- na.omit(data_sub[data_sub$run == 2, ])

## Scale x and y
data_sub3$s_onset_rel    <- scale_m0_sd0.5(data_sub3$onset_rel)
data_sub3$s_noveltyScore <- scale(data_sub3$noveltyScore)

## Fit brms model
m_noveltyScore_run2 <- brm(bf(s_noveltyScore  ~ s_onset_rel  + (s_onset_rel | subject)), 
                           data = data_sub3,
                           family = student(), 
                           prior = prior_noveltyScore, 
                           iter = iter,
                           chains = chains,
                           cores = cores,
                           seed = 1125,
                           control = list(adapt_delta = .95),
                           silent = 0,
                           refresh = 1,
                           backend = "cmdstanr")

summary(m_noveltyScore_run2)

# Save the model
save(m_noveltyScore_run1, m_noveltyScore_run2, file = "intermediate_data/SpaNov_contPM_smo6_m_noveltyScore.Rdata")

