# Script that calculates the slope for the novelty gradient of each subject
# Date:  09/07/2025
# Author: Joern Alexander Quent
# /* 
# ----------------------------- Input for the analysis ---------------------------
# */
# Example path "/home1/Jaquent/Datasets/VP00157_HCP/sessions/R0015_7TSESS1/hcp/R0015_7TSESS1/MNINonLinear/Results/tfMRI_OLMe/tfMRI_OLMe_SpaNov_gradient_6lvl_cue-delay_smo4_MSMAll_hp200_s4_level2_MSMAll_hp0_clean.feat/GrayordinatesStats/cope1.feat"
cope_example_path <- "/home1/Jaquent/Datasets/VP00157_HCP/sessions/XXXX_7TSESS1/hcp/XXXX_7TSESS1/MNINonLinear/Results/tfMRI_OLMe/tfMRI_OLMe_SpaNov_gradient_6lvl_cue-delay_smo4_MSMAll_hp200_s4_level2_MSMAll_hp0_clean.feat/GrayordinatesStats/cope"
GS_example_path   <- "/home1/Jaquent/Datasets/VP00157_HCP/sessions/XXXX_7TSESS1/hcp/XXXX_7TSESS1/MNINonLinear/Results/tfMRI_OLMe_YYY_AP/tfMRI_OLMe_YYY_AP_Atlas_MSMAll_hp0_clean_GScalc.dscalar.nii"
tSNR_example_path <- "/home1/Jaquent/Datasets/VP00157_HCP/sessions/XXXX_7TSESS1/hcp/XXXX_7TSESS1/MNINonLinear/Results/tfMRI_OLMe_YYY_AP/tfMRI_OLMe_YYY_AP_Atlas_MSMAll_hp0_clean_TSNRcalc.dscalar.nii"

# Patterns
subject_pattern <- "XXXX"
run_pattern     <- "YYY"

# Contrasts
copeNums  <- 8:13
conN      <- length(copeNums)
novLabels <- paste0('lvl', 1:conN)

# /* 
# ----------------------------- Paths & other parameters fixed ---------------------------
# */
# Paths
wb_path          <- '/home1/Jaquent/Toolboxes/workbench/bin_rh_linux64/'
sessions_folder  <- '/home1/Jaquent/Datasets/VP00157_HCP/sessions/'
path2subjectFile <- '/home1/Jaquent/research_projects/OLM_project/analsyis/imaging/scripts/SpaNov/subjects2analyse_SpaNov.txt'
path2output      <- '/home1/Jaquent/Datasets/VP00157_HCP/analysis/extracted_values'
path2MNIcoord    <- '/home1/Jaquent/Atlases/cifti_subcortical_MNI152_coordinates.csv'
surfLeft         <- "/home1/Jaquent/Atlases/MMP1.0_210V_Parcellation/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii"
surfRight        <- "/home1/Jaquent/Atlases/MMP1.0_210V_Parcellation/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii"
parcellationFile <- "/home1/Jaquent/Atlases/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors_with_Atlas_ROIs2.32k_fs_LR.dlabel.nii"

# /* 
# ----------------------------- Load libraries etc. ---------------------------
# */
library(stringr)
library(assortedRFunctions)
library(viridis)
library(plyr)

# Load ciftiTools and set workbench paths
possible_wb_paths <- c("/usr/bin/wb_command", "/home1/Jaquent/Toolboxes/workbench/bin_rh_linux64/")
load_ciftiTools(possible_wb_paths)

# /* 
# ----------------------------- Files for ciftiTools etc. ---------------------------
# */
# Load the coordinates
MNI_coord <- read.csv(path2MNIcoord)
load("/home1/Jaquent/Atlases/projected_HC.RData")

# /* 
# ----------------------------- Subject for which to extract values ---------------------------
# */
# Read subject file
subjectFile <- readLines(path2subjectFile)
subjects <- str_split(subjectFile, pattern = ",")[[1]] 
n <- length(subjects)


# /* 
# ----------------------------- Functions ---------------------------
# */
# Function to create gradient xii variables for min and max
create_gradient_xiis <- function(conN, folderName, copeNums, fileName, keyLabels, conColours){
  # Novelty gradient folders
  grad_folders <- paste0(folderName, copeNums)
  
  # Select the correct GLM
  xii_list   <- list()
  
  # Create inverted levels for minimum
  invertedLevels <- as.integer(conN:1)
  
  # Loop through all files
  for(i in 1:length(grad_folders)){
    # Get the current folder 
    currentFolder <- grad_folders[i]
    
    # Load CIFTI file
    cifti_fname   <- paste0(currentFolder, fileName)
    xii_list[[i]] <- ciftiTools::read_cifti(cifti_fname, brainstructures = "all", 
                                            surfL_fname = surfLeft, 
                                            surfR_fname = surfRight)
  }
  
  # Create gradient map by using different label map
  cifti_fname    <- parcellationFile
  grad_max_xii <- ciftiTools::read_cifti(cifti_fname, brainstructures = "all",
                                         surfL_fname = surfLeft,
                                         surfR_fname = surfRight)
  
  # Replace everything with a zero
  grad_max_xii$data$cortex_left  <- matrix(as.integer(rep(0, length(grad_max_xii$data$cortex_left))), ncol = 1)
  grad_max_xii$data$cortex_right <- matrix(as.integer(rep(0, length(grad_max_xii$data$cortex_right))), ncol = 1)
  grad_max_xii$data$subcort      <- matrix(as.integer(rep(0, length(grad_max_xii$data$subcort))), ncol = 1)
  
  # At the same also create the minimum gradient
  grad_min_xii <- grad_max_xii
  
  # Loop through all vertices and voxels and assign label based on zstat
  ## cortex_left
  for(j in 1:length(grad_max_xii$data$cortex_left)){
    # Initialise values
    values <- c()
    
    # Loop through all levels
    for(i in 1:conN){
      values[i] <- xii_list[[i]]$data$cortex_left[j]
    }
    
    # Change the labels in the xii vars
    grad_max_xii$data$cortex_left[j] <- as.integer(which(values == max(values))[1])
    grad_min_xii$data$cortex_left[j] <- invertedLevels[as.integer(which(values == min(values))[1])]
    
    # Progress
    #cat(paste0("\r", j, "/", length(grad_max_xii$data$cortex_left)))
  }
  
  
  ## cortex_right
  for(j in 1:length(grad_max_xii$data$cortex_right)){
    # Initialise values
    values <- c()
    
    # Loop through all levels
    for(i in 1:conN){
      values[i] <- xii_list[[i]]$data$cortex_right[j]
    }
    
    # Change the labels in the xii vars
    grad_max_xii$data$cortex_right[j] <- as.integer(which(values == max(values))[1])
    grad_min_xii$data$cortex_right[j] <- invertedLevels[as.integer(which(values == min(values))[1])]
    
    # Progress
    #cat(paste0("\r", j, "/", length(grad_max_xii$data$cortex_right)))
  }
  
  ## subcort
  for(j in 1:length(grad_max_xii$data$subcort)){
    # Initialise values
    values <- c()
    
    # Loop through all levels
    for(i in 1:conN){
      values[i] <- xii_list[[i]]$data$subcort[j]
    }
    
    # Change the labels in the xii vars
    grad_max_xii$data$subcort[j] <- as.integer(which(values == max(values))[1])
    grad_min_xii$data$subcort[j] <- invertedLevels[as.integer(which(values == min(values))[1])]
  }
  
  # New colours
  conColours <- viridis(n = conN, option = "C")
  
  # Get the old labels
  old_key_colours <- grad_min_xii$meta$cifti$labels$`vertex areas`
  
  # Create new colours
  new_key_colours <- old_key_colours[-(2:nrow(old_key_colours)),]
  
  for(i in 1:conN){
    # Set the colours
    RGB_Col <- col2rgb(conColours[i])/255
    temp_key_colour <- data.frame(Key = i,
                                  Red = RGB_Col[1],
                                  Green = RGB_Col[2],
                                  Blue = RGB_Col[3],
                                  Alpha = 1)
    
    new_key_colours <- rbind(new_key_colours, temp_key_colour)
  }
  
  # Add key labels
  row.names(new_key_colours) <- c("???", keyLabels)
  
  # Add back to xii
  grad_min_xii$meta$cifti$labels$`vertex areas` <- new_key_colours
  grad_max_xii$meta$cifti$labels$`vertex areas` <- new_key_colours
  
  # Save in list
  return(list(grad_min_xii = grad_min_xii, grad_max_xii = grad_max_xii))
}

get_HC_data_from_xii <- function(currentxii){
  # Extract hippocampal data and combine with MNI coordinates
  # Left
  HC_data_left     <- MNI_coord[str_starts(MNI_coord$region, pattern = "Hippocampus-L"), ]
  index_bool       <- str_starts(currentxii$grad_min_xii$meta$subcort$labels, pattern = "Hippocampus-L")
  HC_data_left$min <- currentxii$grad_min_xii$data$subcort[index_bool, ]
  HC_data_left$max <- currentxii$grad_max_xii$data$subcort[index_bool, ]
  HC_data_left$Hemisphere <- "left"
  
  # Add new column/variable for position
  HC_data_left$position <- NA
  
  # Add the projected position
  projected_HC_L <- projected_HC[projected_HC$region == "Hippocampus-L", ]
  
  for(i in 1:nrow(projected_HC_L)){
    # Get the values for the projection
    y <- projected_HC_L$y[i]
    z <- projected_HC_L$z[i]
    pos <- projected_HC_L$position[i]
    
    # Add position values
    HC_data_left[HC_data_left$y == y & HC_data_left$z == z, "position"] <- pos
  }
  
  # Right
  HC_data_right     <- MNI_coord[str_starts(MNI_coord$region, pattern = "Hippocampus-R"), ]
  index_bool        <- str_starts(currentxii$grad_min_xii$meta$subcort$labels, pattern = "Hippocampus-R")
  HC_data_right$min <- currentxii$grad_min_xii$data$subcort[index_bool, ]
  HC_data_right$max <- currentxii$grad_max_xii$data$subcort[index_bool, ]
  HC_data_right$Hemisphere <- "right"
  
  # Add new column/variable for position
  HC_data_right$position <- NA
  
  # Add the projected position
  projected_HC_R <- projected_HC[projected_HC$region == "Hippocampus-R", ]
  
  for(i in 1:nrow(projected_HC_R)){
    # Get the values for the projection
    y <- projected_HC_R$y[i]
    z <- projected_HC_R$z[i]
    pos <- projected_HC_R$position[i]
    
    # Add position values
    HC_data_right[HC_data_right$y == y & HC_data_right$z == z, "position"] <- pos
  }
  
  # Concatenate to one data frame
  HC_data <- rbind(HC_data_left, HC_data_right)
  
  # Return
  return(HC_data)
}

get_GS_and_tSNR_values <- function(subject){
  # Load GS values
  RUN1_cifti_file <- str_replace_all(GS_example_path, pattern = subject_pattern, replacement = subject)
  RUN2_cifti_file <- str_replace_all(GS_example_path, pattern = subject_pattern, replacement = subject)
  RUN1_cifti_file <- str_replace_all(RUN1_cifti_file, pattern = run_pattern, replacement = "RUN1")
  RUN2_cifti_file <- str_replace_all(RUN2_cifti_file, pattern = run_pattern, replacement = "RUN1")
  RUN1_xii        <- read_cifti(RUN1_cifti_file, brainstructures = "all")
  RUN2_xii        <- read_cifti(RUN2_cifti_file, brainstructures = "all")
  
  ## Calculate average across both runs
  avg_GS <- (RUN1_xii + RUN2_xii)/2
  
  # Load tSNR values
  RUN1_cifti_file <- str_replace_all(tSNR_example_path, pattern = subject_pattern, replacement = subject)
  RUN2_cifti_file <- str_replace_all(tSNR_example_path, pattern = subject_pattern, replacement = subject)
  RUN1_cifti_file <- str_replace_all(RUN1_cifti_file, pattern = run_pattern, replacement = "RUN1")
  RUN2_cifti_file <- str_replace_all(RUN2_cifti_file, pattern = run_pattern, replacement = "RUN1")
  RUN1_xii        <- read_cifti(RUN1_cifti_file, brainstructures = "all")
  RUN2_xii        <- read_cifti(RUN2_cifti_file, brainstructures = "all")
  
  ## Calculate average across both runs
  avg_tSNR <- (RUN1_xii + RUN2_xii)/2
  
  # Get bool index for left & right HC
  HC_index_L <- avg_tSNR$meta$subcort$labels == "Hippocampus-L"
  HC_index_R <- avg_tSNR$meta$subcort$labels == "Hippocampus-R"
  
  # Extract hippocampus values
  GS_values    <- c(avg_GS$data$subcort[HC_index_L], avg_GS$data$subcort[HC_index_R])
  tSNR_values  <- c(avg_tSNR$data$subcort[HC_index_L], avg_tSNR$data$subcort[HC_index_R])
  
  # Return as data frame
  return(data.frame(GS = GS_values, tSNR = tSNR_values))
}

# /* 
# ----------------------------- Loop through all subjects ---------------------------
# */
# Prepare data frame
subjSlopes <- as.data.frame(matrix(nrow = 0, ncol = 7))

# Start time for progress
startTime <- Sys.time()

# Loop through everyone
for(subj in 1:length(subjects)){
  # Combine subject with suffix
  subject_string <- subjects[subj]
  
  # Create folder name by combining all parts
  folderName <-str_replace_all(cope_example_path, pattern = subject_pattern, replacement = subject_string)

  # Create the current gradient xifti for this subject
  currentxii <- create_gradient_xiis(conN, folderName, copeNums, ".feat/cope1.dtseries.nii", novLabels)

  # Get the hippocampus data
  HC_data <- get_HC_data_from_xii(currentxii)
  
  # Get GS & tSNR for hippocampus and add to HC_data
  co_variates <- get_GS_and_tSNR_values(subject_string)
  HC_data <- cbind(HC_data, co_variates)

  # Average across voxels
  HC_data_agg1_y <- ddply(HC_data, c("y", "Hemisphere"), summarise, min = mean(min), max = mean(max), 
                          tSNR = mean(tSNR), GS = mean(GS))
  HC_data_agg1_pos <- ddply(HC_data, c("position", "Hemisphere"), summarise, n = length(min), min = mean(min), max = mean(max), 
                            tSNR = mean(tSNR), GS = mean(GS))

  # Fit models
  sum1 <- as.data.frame(summary(lm(min ~ y * Hemisphere + GS + tSNR, data = HC_data_agg1_y))$coefficients)
  sum2 <- as.data.frame(summary(lm(min ~ position * Hemisphere + GS + tSNR, data = HC_data_agg1_pos))$coefficients)
  
  # Subset to left & right 
  HC_data_agg1_y_L   <- HC_data_agg1_y[HC_data_agg1_y$Hemisphere == "left", ]
  HC_data_agg1_y_R   <- HC_data_agg1_y[HC_data_agg1_y$Hemisphere == "right", ]
  HC_data_agg1_pos_L <- HC_data_agg1_pos[HC_data_agg1_pos$Hemisphere == "left", ]
  HC_data_agg1_pos_R <- HC_data_agg1_pos[HC_data_agg1_pos$Hemisphere == "right", ]
  
  # Fit models separately for each hemisphere
  sum3 <- as.data.frame(summary(lm(min ~ y + GS + tSNR, data = HC_data_agg1_y_L))$coefficients)
  sum4 <- as.data.frame(summary(lm(min ~ y + GS + tSNR, data = HC_data_agg1_y_R))$coefficients)
  sum5 <- as.data.frame(summary(lm(min ~ position + GS + tSNR, data = HC_data_agg1_pos_L))$coefficients)
  sum6 <- as.data.frame(summary(lm(min ~ position + GS + tSNR, data = HC_data_agg1_pos_R))$coefficients)
  
  # Add extra information
  sum1$subject <- subjects[subj]
  sum1$type    <- "MNIy"
  sum1$effect  <- row.names(sum1)
  
  sum2$subject <- subjects[subj]
  sum2$type    <- "position"
  sum2$effect  <- row.names(sum2)
  
  sum3$subject <- subjects[subj]
  sum3$type    <- "MNIy_L"
  sum3$effect  <- row.names(sum3)
  
  sum4$subject <- subjects[subj]
  sum4$type    <- "MNIy_R"
  sum4$effect  <- row.names(sum4)
  
  sum5$subject <- subjects[subj]
  sum5$type    <- "position_L"
  sum5$effect  <- row.names(sum5)
  
  sum6$subject <- subjects[subj]
  sum6$type    <- "position_R"
  sum6$effect  <- row.names(sum6)

  # Add to results
  subjSlopes <- rbind(subjSlopes, sum1, sum2, sum3, sum4, sum5, sum6)

  # Progress
  progressDisplay(subj, n, startTime)
}


# Remove row names
row.names(subjSlopes) <- NULL

# /* 
# ----------------------------- Save data frame to disk---------------------------
# */
cat('\n')
outputFileName <- datedFileNam(fileName = "SpavNov_gradient_subject-slopes", fileEnding = '.RData')
dir.create(path2output, showWarnings = FALSE, recursive = TRUE)
save(subjSlopes, file = paste(path2output, outputFileName, sep = '/'))


