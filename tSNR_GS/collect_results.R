# Script to collect the GS & tSNR maps
# Date:  06/07/2025
# Author: Joern Alexander Quent
# /* 
# ----------------------------- Libraries ---------------------------
# */
# Libraries
library(stringr)
library(assortedRFunctions)

# Load ciftiTools and set workbench paths
possible_wb_paths <- c("/usr/bin/wb_command", "/home1/Jaquent/Toolboxes/workbench/bin_rh_linux64/")
load_ciftiTools(possible_wb_paths)

# /* 
# ----------------------------- Load the paths ---------------------------
# */
# Load paths
GS_paths   <- read.csv("GS_paths.txt", header = FALSE)$V1
tSNR_paths <- read.csv("tSNR_paths.txt", header = FALSE)$V1

# Subset to OLMe 
GS_paths   <- GS_paths[str_detect(GS_paths, pattern = "OLMe")]
tSNR_paths <- tSNR_paths[str_detect(tSNR_paths, pattern = "OLMe")]

# /* 
# ----------------------------- Loop through all images ---------------------------
# */
GS_list      <- list()
tSNR_list    <- list()
GS_subject   <- c()
GS_run       <- c()
tSNR_subject <- c()
tSNR_run     <- c()

for(i in 1:length(GS_paths)){
  # Print console
  cat("\nProcessing", i, "out of", length(GS_paths))
  
  # GS
  ## Get path
  path <- GS_paths[i]
  
  ## Parse information from path 
  GS_subject[i] <- str_split_1(str_split_i(path, pattern = "/", 7), pattern = "_")[1]
  GS_run[i]     <- str_split_1(str_split_i(path, pattern = "/", 13), pattern = "_")[3]
  
  ## Load cifti
  GS_list[[i]] <- read_cifti(path, brainstructures = "all")
  
  #tSNR
  ## Get path
  path <- tSNR_paths[i]
  
  ## Parse information from path 
  tSNR_subject[i] <- str_split_1(str_split_i(path, pattern = "/", 7), pattern = "_")[1]
  tSNR_run[i]     <- str_split_1(str_split_i(path, pattern = "/", 13), pattern = "_")[3]
  
  ## Load cifti
  tSNR_list[[i]] <- read_cifti(path, brainstructures = "all")
}

# Create data frame 
GS_tSNR_df <- data.frame(GS_subject, GS_run, tSNR_subject, tSNR_run)

# Save variables
save(GS_tSNR_df, GS_list, tSNR_list, file = "tSNR_GS_maps.RData")