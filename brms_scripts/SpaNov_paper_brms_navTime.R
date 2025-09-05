# Script to run the brms model for navTime analysis of the SpaNov paper
# Date:  15/07/2025
# Author: Joern Alexander Quent
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
library(data.table)
library(stringr)
library(assortedRFunctions)

# /* 
# ----------------------------- Load data and prepare ---------------------------
# */
# Specify paths where the data is saved
path2data <- "data/ignore_fMRI_version1/"

# Load data
load(paste0(path2data, "combined_data/OLM_7T_all_data.RData"))

# Load the look-up table that contains information of R-numbers which are retracted 
lookupTable  <- read.csv(paste0(path2data, "lookUpTable.csv"))

# Select the subjects included in this analysis
## Load the subjects that are included in this analysis
subjectFile <- readLines(paste0(path2data, "SpaNov_subject2analyse.txt"))
subjIDs_R   <- str_split(subjectFile, pattern = ",")[[1]] 
subjIDs     <- lookupTable$anonKey[lookupTable$Rnum %in% subjIDs_R]

OLM_7T_trial_results <- OLM_7T_trial_results[OLM_7T_trial_results$subject %in% subjIDs, ]

# Subset to encoding trials
## Use unique name for data frame because of marginaleffects issues
m_navTime_data <- OLM_7T_trial_results[OLM_7T_trial_results$trialType == "encoding", ]

# Calculate distance between start and object
x1 <- m_navTime_data$start_x
z1 <- m_navTime_data$start_z
x2 <- m_navTime_data$object_x
z2 <- m_navTime_data$object_z
m_navTime_data$start2obj <- euclideanDistance3D(x1, 1, z1, x2, 1, z2)

# Scale values
m_navTime_data$s_start2obj <- scale_m0_sd0.5(m_navTime_data$start2obj)
m_navTime_data$s_timesObjectPresented <- scale_m0_sd0.5(m_navTime_data$timesObjectPresented)

# /* 
# ----------------------------- s_navTime ~ s_start2obj  ---------------------------
# */
# Set brms priors
priors  <- c(prior(normal(0, 1), class = "Intercept"),
             prior(normal(0, 1), class = "b"))

# Fit brms model
m_navTime1 <- brm(navTime ~ s_start2obj + (s_start2obj | subject), 
               data = m_navTime_data,
               family = lognormal(), 
               prior = priors, 
               iter = iter,
               chains = chains, 
               cores = cores, 
               seed = 1421,
               #control = list(adapt_delta = .95),
               silent = 0,
               refresh = 1,
               save_pars = save_pars(all = TRUE),
               backend = "cmdstanr")

# Show results
summary(m_navTime1)

# /* 
# ----------------------------- s_navTime ~ s_timesObjectPresented  ---------------------------
# */
# Set brms priors
priors  <- c(prior(normal(0, 1), class = "Intercept"),
             prior(normal(0, 1), class = "b"))

# Fit brms model
m_navTime2 <- brm(navTime ~ s_timesObjectPresented + (s_timesObjectPresented | subject), 
                  data = m_navTime_data,
                  family = lognormal(), 
                  prior = priors, 
                  iter = iter,
                  chains = chains, 
                  cores = cores, 
                  seed = 1421,
                  #control = list(adapt_delta = .95),
                  silent = 0,
                  refresh = 1,
                  save_pars = save_pars(all = TRUE),
                  backend = "cmdstanr")

# Show results
summary(m_navTime2)

# /* 
# ----------------------------- navTime ~ s_start2obj + s_timesObjectPresented ---------------------------
# */
# Set brms priors
priors  <- c(prior(normal(0, 1), class = "Intercept"),
             prior(normal(0, 1), class = "b"))

# Fit brms model
m_navTime3 <- brm(navTime ~ s_start2obj + s_timesObjectPresented + (s_start2obj + s_timesObjectPresented | subject), 
                  data = m_navTime_data,
                  family = lognormal(), 
                  prior = priors, 
                  iter = iter,
                  chains = chains, 
                  cores = cores, 
                  seed = 1421,
                  #control = list(adapt_delta = .95),
                  silent = 0,
                  refresh = 1,
                  save_pars = save_pars(all = TRUE),
                  backend = "cmdstanr")

# Show results
summary(m_navTime3)


# /* 
# ----------------------------- Save results ---------------------------
# */
# Save the model
save(m_navTime1, m_navTime2, m_navTime3, file = "fitted_brms_models/SpaNov_m_navTime.Rdata")
