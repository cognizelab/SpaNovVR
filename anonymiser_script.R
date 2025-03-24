# This script serves two functions a) anonymising the subject IDs & b) create homogenise subject IDs across tasks

# Set wd
#setwd("~/GitHub/OLM_project/analysis")
setwd("~/Work/researchProjects/OLM_project/analysis")

# Libs
library(assortedRFunctions)
library(readxl)

# Load data collection spreadsheet
spreadsheet <- read_excel("data/ignore_fMRI_version1/VP00157_PPT_Info_Final.xlsx")
#spreadsheet <- spreadsheet[-1,] # Remove first row because this is only an example

# Create look up table
lookUp <- data.frame(Rnum = spreadsheet$`ID Number`, HCnum = spreadsheet$`P Number`)

# Set seed & anonymise
set.seed(20221114) # This will ensure that only the participants that will be added get a new number
lookUp$anonKey <- anonymise(lookUp$Rnum, fileName = "data/ignore_fMRI_version1/anonKey")

# Also save the look out table
write.csv(lookUp, "data/ignore_fMRI_version1/lookUpTable.csv", row.names = FALSE)