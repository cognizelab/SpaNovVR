# Script to run the functional resting state connectivity analysis between each hippocampal voxel and each cortical vertex
# Note: This script is supposed to be run on the cluster. 
# Date:  18/11/2024
# Author: Joern Alexander Quent
# /* 
# ----------------------------- Preparation ---------------------------
# */
# Libraries
library(foreach)
library(doParallel)
library(stringr)
library(ciftiTools)

# Paths
wb_path <- '/home1/Jaquent/Toolboxes/workbench/bin_rh_linux64/'

# Surface files
surfLeft         <- "../../source_files/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii"
surfRight        <- "../../source_files/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii"

# Load the R-numbers of the subjects that are supposed to be analysed.
sub2analyse_file <-"../SpaNov/subjects2analyse_SpaNov.txt"
sub2analyse      <- readLines(sub2analyse_file)
sub2analyse      <- str_split(sub2analyse, pattern = ",")[[1]]
subjects_num     <- length(sub2analyse)

# Print
cat("\nA total of", length(sub2analyse), "subjects were in the subject list.\n")

# Base path to the resting state data
pattern      <- "XXXXX" # Replace this with the correct subject R-number.
basePath_RSD <- "/home1/Jaquent/Datasets/VP00157_HCP/sessions/XXXXX_3TSESS1/hcp/XXXXX_3TSESS1/MNINonLinear/Results/REST/REST_Atlas_MSMAll_hp0_clean.dtseries.nii"

# Base name of output file
output_fileName <- "data/HC_2_cortex_FRSC_XXXXX.RData"

# Create custom log
logFile <- "log/job_dopar.out"
writeLines(c(""), logFile)

# Write to log
sink(logFile, append=TRUE)
cat("\n This job has", detectCores(), "CPUs.\n")

# /* 
# ----------------------------- Process the resting state data and get the time series ---------------------------
# */
# Start cluster
my.cluster <- parallel::makeCluster(detectCores() - 2, type = "PSOCK")

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Loop through all subjects
HC_2_cortex_timeseries <- foreach(subj = 1:subjects_num, .packages = c('stringr', 'ciftiTools')) %dopar% {
  # Start writing to the log file
  sink(logFile, append=TRUE)
  
  # Get the current R-number
  Rnum <-sub2analyse[subj]
  
  # Set workbench path
  ciftiTools.setOption('wb_path', wb_path)
  
  # Load the RSD from the current subject
  temp_path <- str_replace_all(basePath_RSD, pattern = pattern, replacement = Rnum)
  temp_xii  <- read_cifti(cifti_fname = temp_path, brainstructures = "all")
  
  # Create matrices for hippocampus and cortex
  HC_VoxelNum     <- sum(str_detect(temp_xii$meta$subcort$labels, pattern = "Hippocampus"))
  temp_HC_mat     <- temp_xii$data$subcort[str_detect(temp_xii$meta$subcort$labels, pattern = "Hippocampus"), ]
  temp_cortex_mat <- rbind(temp_xii$data$cortex_left, temp_xii$data$cortex_right)
  # The left cortex vertices should be at the top, right cortex vertices in the middle, and subcortex vertices at the bottom (when present).
  
  # Create empty matrix with zeros because it will include the subcortical values even though they will be ignored. 
  # Rows = all grayordinates, columns = hippocampal voxels
  conn_mat <- matrix(0, nrow = nrow(temp_cortex_mat) + nrow(temp_xii$data$subcort),
                     ncol = HC_VoxelNum)
  
  # Loop through each cortical vertex
  for(i in 1:nrow(temp_cortex_mat)){
    # Get cortical time series
    t1 <- temp_cortex_mat[i, ]
    
    # Loop through each hippocampal voxel
    for(j in 1:HC_VoxelNum){
      # Get hippocampal time series
      t2 <- temp_HC_mat[j, ]
      
      # Add to connectivity matrix
      conn_mat[i, j] <- cor(t1, t2)
    }
    
    # Print progress every 500 vertices
    if(i %% 500 == 0){
      cat("\nProcessing", Rnum, paste0("(", subj, "/", length(sub2analyse), ")"), 
          "|" , as.character(Sys.time()), "|", 
          round((i/nrow(temp_cortex_mat)) * 100, 3), "% completed") 
    }
  }
  
  # Create custom file name
  tempName <- str_replace_all(output_fileName, pattern = pattern, replacement = Rnum)
  
  # Print name
  cat("\nWriting ", tempName)
  
  # Write matrix to file
  save(conn_mat, file = tempName)
}
