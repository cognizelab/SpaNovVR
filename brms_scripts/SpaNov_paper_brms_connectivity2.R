# Script to run the brms model for connectivity analysis of the SpaNov paper
# Date:  07/07/2025
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
library(plyr)

# /* 
# ----------------------------- Load data and prepare ---------------------------
# */
# Load FRSC data
load("data/ignore_fMRI_version1/grad_FRSC/HC_2_cortex_FRSC_SpaNov_gradient_cue-delay.RData")

# Convert list to data frame
## Use unique name for data frame because of marginaleffects issues
m_conn2_data <- rbindlist(HC_2_cortex_SpaNov_gradient)
m_conn2_data <- as.data.frame(m_conn2_data)

## Average across hemipsheres
m_conn2_data <- ddply(m_conn2_data, c("subject", "gradient_level_cortex", "gradient_level_HC"), summarise,
                      connectivity  = mean(connectivity, na.rm = TRUE))


# Create column whether or not a pair is on the diagonal
m_conn2_data$diagonal <- ifelse(m_conn2_data$gradient_level_cortex == m_conn2_data$gradient_level_HC, "diagonal", "off-diagonal")

# Create other factors
m_conn2_data$f_gradient_level_cortex <- factor(m_conn2_data$gradient_level_cortex)
m_conn2_data$f_gradient_level_HC     <- factor(m_conn2_data$gradient_level_HC)

# Scale variables
m_conn2_data$s_connectivity <- drop(scale(m_conn2_data$connectivity))


# /* 
# ----------------------------- connectivity ~ gradient_level_cortex + gradient_level_HC + diagonal ---------------------------
# */
# Set brms priors
priors  <- c(prior(normal(0, 1), class = "Intercept"),
             prior(normal(0, 1), class = "b"))

# Fit brms model
m_conn2 <- brm(s_connectivity ~ f_gradient_level_cortex + f_gradient_level_HC + diagonal + 
                             (f_gradient_level_cortex + f_gradient_level_HC + diagonal | subject), 
               data = m_conn2_data,
               family = gaussian(), 
               prior = priors, 
               iter = iter,
               chains = chains, 
               cores = cores, 
               seed = 31231,
               #control = list(adapt_delta = .95),
               silent = 0,
               refresh = 1,
               backend = "cmdstanr")

# Show results
summary(m_conn2)

# Save the model
save(m_conn2, file = "fitted_brms_models/SpaNov_m_conn2.Rdata")
