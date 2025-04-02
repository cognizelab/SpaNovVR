# Script that analyses the estimated HC 2 cortex connectivity using the spatial novelty gradients.
# Note: This script is supposed to be run on the cluster. 
# Date:  20/11/2024
# Author: Joern Alexander Quent
# /* 
# ----------------------------- Libraries & parameters ---------------------------
# */
# Libraries
library(assortedRFunctions)
library(plyr)
library(stringr)
library(ciftiTools)
library(foreach)
library(doParallel)

# Paths
wb_path <- '/home1/Jaquent/Toolboxes/workbench/bin_rh_linux64/'
ciftiTools.setOption('wb_path', wb_path)

# Job parameters
output_fileName <- "results/HC_2_cortex_FRSC_SpaNov_gradient.RData"
input_fileName  <- "data/HC_2_cortex_FRSC_XXXXX.RData"
job_name        <- "SpaNov_gradient"
pattern         <- "XXXXX" # Replace this with the correct subject R-number.


# Create custom log
logFile <- paste0("log/job_", job_name, ".out")

# Write to log
sink(logFile, append=TRUE)
cat("\n This job has", detectCores(), "CPUs.\n")

# /* 
# ----------------------------- Load and prepare ---------------------------
# */
# Load the R-numbers of the subjects that are supposed to be analysed.
sub2analyse_file <-"../SpaNov/subjects2analyse_SpaNov.txt"
sub2analyse      <- readLines(sub2analyse_file)
sub2analyse      <- str_split(sub2analyse, pattern = ",")[[1]]
subjects_num     <- length(sub2analyse)

# Load the CIFTI files and prepare for the analysis
## Load the CIFTI files
PMC_gradient         <- read_cifti("gradient_maps/SpaNov_gradient_PMC_min.dlabel.nii", brainstructures = "all")
whole_brain_gradient <- read_cifti("gradient_maps/SpaNov_gradient_wholebrain_min.dlabel.nii", brainstructures = "all")

## Prepare vector with labels
gradient_level_cortex <- get_all_points_from_xifti(PMC_gradient) # left, right & subcortical
hemisphere_cortex     <- c(rep("left", length(PMC_gradient$data$cortex_left)),
                    rep("right", length(PMC_gradient$data$cortex_right)),
                    rep("subcort", length(PMC_gradient$data$subcort)))
gradient_level_HC <- whole_brain_gradient$data$subcort[str_detect(whole_brain_gradient$meta$subcort$labels, pattern = "Hippocampus")]
hemisphere_HC     <- whole_brain_gradient$meta$subcort$labels[str_detect(whole_brain_gradient$meta$subcort$labels, pattern = "Hippocampus")]
hemisphere_HC     <- ifelse(as.character(hemisphere_HC) == "Hippocampus-R", "right", "left")

# Prepare data frame to save data by combining all levels
template_DF <- expand.grid(hemisphere_cortex = c("left", "right"),
                           hemisphere_HC = c("left", "right"),
                           gradient_level_cortex = 1:6,
                           gradient_level_HC = 1:6)
template_DF$connectivity <- NA

# /* 
# ----------------------------- Load subject connectivity matrices and aggregate based on the novelty levels ---------------------------
# */
# Start cluster
my.cluster <- parallel::makeCluster(detectCores(), type = "PSOCK")

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Loop through all subjects
HC_2_cortex_SpaNov_gradient <- foreach(subj = 1:subjects_num, .packages = c("stringr", "assortedRFunctions")) %dopar% {
  # Start writing to the log file
  sink(logFile, append=TRUE)
  
  # Get the current R-number
  Rnum <-sub2analyse[subj]
  
  # Add R-number to template_DF
  template_DF$subject <- Rnum
  
  # Get path to data
  temp_path <- str_replace_all(input_fileName, pattern = pattern, replacement = Rnum)
  
  # Load participant matrix
  load(temp_path)
  
  # Apply fisher z transformation to the correlation matrix
  conn_mat <- z_transform_fisher(conn_mat)
  
  # Calculate average connectivity between PMC and hippocampus for each novelty level
  # by looping through template_DF
  for(i in 1:nrow(template_DF)){
    # Get indices for the rows & columns
    rowID <- which(gradient_level_cortex == template_DF$gradient_level_cortex[i] & hemisphere_cortex == template_DF$hemisphere_cortex[i])
    colID <- which(gradient_level_HC == template_DF$gradient_level_HC[i] & hemisphere_HC == template_DF$hemisphere_HC[i])
    
    # Calculate mean connectivity
    template_DF$connectivity[i] <- mean(conn_mat[rowID, colID])
  }
  
  cat("\nProcessing", Rnum, paste0("(", sprintf("%03d", subj), "/", sprintf("%03d", length(sub2analyse)), ")"), 
      "|" , as.character(Sys.time()), "\n") 
  
  # Return value
  template_DF
}

# /* 
# ----------------------------- Load subject connectivity matrices and aggregate based on the novelty levels ---------------------------
# */
save(HC_2_cortex_SpaNov_gradient, file = output_fileName)
