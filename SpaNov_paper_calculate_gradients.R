# Script to calculate the gradients for the SpaNov paper to clean up the main markdown file
# Date:  25/11/2024
# Author: Joern Alexander Quent
# /* 
# ----------------------------- Libraries ---------------------------
# */
# Libraries
library(ciftiTools)

# Use correct locations and other settings based on computer
if(Sys.info()[4] == "DESKTOP-335I26I"){
  # Work laptop (Windows)
  ## Setting paths to workbench installation
  ciftiTools.setOption("wb_path", "C:/Program Files/workbench-windows64-v1.5.0/workbench/bin_windows64")
  
  ## Path to the imaging data
  path2imaging_results2 <- "D:/Seafile/imaging_results"
  
} else if(Sys.info()[4] == 'DESKTOP-91CQCSQ') {
  # Work desktop (Windows)
  ## Setting paths to workbench installation
  ciftiTools.setOption("wb_path", "D:/workbench/bin_windows64")
  
  ## Path to the imaging data
  path2imaging_results2 <- "D:/imaging_results"
} else if(Sys.info()[4] == 'alex-Zenbook-UX3404VA-UX3404VA') {
  # Work laptop (Linux)
  ## Setting paths to workbench installation
  ciftiTools.setOption("wb_path", "/usr/bin/wb_command")
  
  ## Path to the imaging data
  path2imaging_results2 <- "/media/alex/shared/Seafile/imaging_results"
  
} else if(Sys.info()[4] == "greengoblin"){
  # Main desktop PC (Linux)
  ciftiTools.setOption("wb_path", "/usr/bin/wb_command") 
  path2imaging_results2 <- "/media/alex/work/Seafile/imaging_results" 
} else if(Sys.info()[4] == "GREEN-GOBLIN-WI"){
  # Main desktop PC (Linux)
  ciftiTools.setOption("wb_path", "C:/Program Files/workbench/bin_windows64")
  path2imaging_results2 <- "E:/Seafile/imaging_results" 
} else {
  # Personal laptop (Windows)
  ## Setting paths to workbench installation
  ciftiTools.setOption("wb_path", "D:/Program Files/workbench/bin_windows64")
  
  ## Path to the imaging data
  path2imaging_results2 <- "D:/OLM/imaging_results"
}

# /* 
# ----------------------------- Stuff get the CFITI files ---------------------------
# */
# Loading the support CIFTI files
## Place where to find some of the CIFTI files and parcellations
CIFTI_locations <- "data/ignore_fMRI_version1/sourceFiles/"

## CIFTI files
parcellationFile <- "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors_with_Atlas_ROIs2.32k_fs_LR.dlabel.nii"
CAB_NP           <- "CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii"
surfLeft         <- "S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii"
surfRight        <- "S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii"

## Combine CIFTI_locations with file name
parcellationFile <- paste0(CIFTI_locations, parcellationFile)
CAB_NP           <- paste0(CIFTI_locations, CAB_NP)
surfLeft         <- paste0(CIFTI_locations, surfLeft)
surfRight        <- paste0(CIFTI_locations, surfRight)

## Loading CIFTIw via ciftiTools as xiis
### Get MMP parcellation
MMP_xii <- ciftiTools::read_cifti(parcellationFile,
                                  brainstructures = "all", 
                                  surfL_fname = surfLeft, 
                                  surfR_fname = surfRight)

### Load Yeo 7 parcellation
Yeo7_xii <- ciftiTools::read_cifti("other_stuff/Yeo7.dlabel.nii",
                                   surfL_fname = surfLeft, 
                                   surfR_fname = surfRight)

## Load other stuff
### Load the parcel names for MMP
parcel_names <- read.csv("data/ignore_fMRI_version1/extracted_values/Parcellations/MP1.0_210V_parcel_names.csv", header = FALSE)

### Load the extracted MNI coordinates
MNI_coord <- read.csv("other_stuff/cifti_subcortical_MNI152_coordinates.csv")

### Load the hippocampal projection values
load("other_stuff/projected_HC.RData")

# /* 
# ----------------------------- Functions used to run the calculations ---------------------------
# */
# Function to create gradient xii variables for min and max
create_gradient_xiis <- function(conN, folderName, copeNums, fileName, keyLabels, mapType = "zstat", conColours){
  # Novelty gradient folders
  grad_folders <- paste0("cope", copeNums, '.feat')
  
  # Select the correct GLM
  xii_list   <- list()
  
  # Create inverted levels for minimum
  invertedLevels <- as.integer(conN:1)
  
  # Loop through all files
  for(i in 1:length(grad_folders)){
    # Get the current folder 
    currentFolder <- paste(path2imaging_results2, folderName, grad_folders[i], sep = "/")
    
    # Load CIFTI file
    cifti_fname   <- paste(currentFolder, fileName, sep = "/")
    if(mapType == 'zstat'){
      xii_list[[i]] <- ciftiTools::read_cifti(cifti_fname, brainstructures = "all", 
                                              surfL_fname = surfLeft, 
                                              surfR_fname = surfRight)
    } else if(mapType == 'beta'){
      xii_beta <- ciftiTools::read_cifti(cifti_fname, brainstructures = "all", 
                                         surfL_fname = surfLeft, 
                                         surfR_fname = surfRight)
      
      # Average over people
      xii_list[[i]] <- apply_xifti(xii_beta, margin = 1, mean)
    }
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


# Taken from https://stackoverflow.com/questions/15882323/r-robust-fitting-of-data-points-to-a-gaussian-function
# Peer et al. (2021) used a different function in matlab (gauss1: Y = a1*exp(-((x-b1)/c1)^2)) but I don't think this matters
fit_gaussian_curve <- function(x, y, mu, sig, scaling){
  f <- function(p){
    d <- p[3]*dnorm(x, mean = p[1],sd = p[2])
    sum((d-y)^2)
  }
  
  optim(c(mu, sig, scaling), f)
}


create_gaussian_gradient_xiis <- function(conN, folderName, copeNums, fileName, keyLabels, mapType = "zstat", conColours){
  # Novelty gradient folders
  grad_folders <- paste0("cope", copeNums, '.feat')
  
  # Select the correct GLM
  xii_list   <- list()
  
  # Create inverted levels for minimum
  invertedLevels <- as.integer(conN:1)
  
  # Loop through all files
  for(i in 1:length(grad_folders)){
    # Get the current folder 
    currentFolder <- paste(path2imaging_results2, folderName, grad_folders[i], sep = "/")
    
    # Load CIFTI file
    cifti_fname   <- paste(currentFolder, fileName, sep = "/")
    if(mapType == 'zstat'){
      xii_list[[i]] <- ciftiTools::read_cifti(cifti_fname, brainstructures = "all", 
                                              surfL_fname = surfLeft, 
                                              surfR_fname = surfRight)
    } else if(mapType == 'beta'){
      xii_beta <- ciftiTools::read_cifti(cifti_fname, brainstructures = "all", 
                                         surfL_fname = surfLeft, 
                                         surfR_fname = surfRight)
      
      # Average over people
      xii_list[[i]] <- apply_xifti(xii_beta, margin = 1, mean)
    }
  }
  
  # Create gradient map by using different label map
  cifti_fname <- parcellationFile
  grad_xii  <- ciftiTools::read_cifti(cifti_fname, brainstructures = "all", 
                                      surfL_fname = surfLeft, 
                                      surfR_fname = surfRight)
  
  # Replace everything with a zero
  grad_xii$data$cortex_left  <- matrix(as.integer(rep(0, length(grad_xii$data$cortex_left))), ncol = 1)
  grad_xii$data$cortex_right <- matrix(as.integer(rep(0, length(grad_xii$data$cortex_right))), ncol = 1)
  grad_xii$data$subcort      <- matrix(as.integer(rep(0, length(grad_xii$data$subcort))), ncol = 1)
  
  # Create scalar xii to save other information from the fitting process
  sigma_xii   <- GLM1_zMap_xii
  
  ## Replace everything with a zero
  sigma_xii$data$cortex_left  <- matrix(rep(0, length(grad_xii$data$cortex_left)), ncol = 1)
  sigma_xii$data$cortex_right <- matrix(rep(0, length(grad_xii$data$cortex_right)), ncol = 1)
  sigma_xii$data$subcort      <- matrix(rep(0, length(grad_xii$data$subcort)), ncol = 1)
  
  ## Use sigma_xii ton create the other
  scaling_xii     <- sigma_xii
  value_xii       <- sigma_xii
  
  # Fixed parameters for the Gaussian fit
  x       <- 1:conN
  mu      <- mean(x)
  sig     <- 3.50114 # Value based the median of run with sd(x)
  scaling <- 31.03736 # Value based the median of run with 10
  
  # Loop through all vertices and voxels and assign label based on zstat
  ## cortex_left
  for(j in 1:length(grad_xii$data$cortex_left)){
    # Initialise values
    values <- c()
    
    # Loop through all levels
    for(i in 1:conN){
      values[i] <- xii_list[[i]]$data$cortex_left[j]
    }
    
    # Make all values positive
    values <- values - min(values)
    
    # Fit Gaussian function
    fit <- fit_gaussian_curve(x, values, mu, sig, scaling)
    
    # Clamp peaks falling outside the range to the range
    fit_mu <- fit$par[1]
    if(fit_mu < 1){
      fit_mu <- 1
    } else if(fit_mu > conN){
      fit_mu <- conN
    }
    
    # Round values to discretise 
    fit_mu <- round(fit_mu)
    
    # Add back to the xii
    grad_xii$data$cortex_left[j]    <- fit_mu
    sigma_xii$data$cortex_left[j]   <- fit$par[2]
    scaling_xii$data$cortex_left[j] <- fit$par[3]
    value_xii$data$cortex_left[j]   <- fit$value
  }
  
  
  ## cortex_right
  for(j in 1:length(grad_xii$data$cortex_right)){
    # Initialise values
    values <- c()
    
    # Loop through all levels
    for(i in 1:conN){
      values[i] <- xii_list[[i]]$data$cortex_right[j]
    }
    
    # Make all values positive
    values <- values - min(values)
    
    # Fit Gaussian function
    fit <- fit_gaussian_curve(x, values, mu, sig, scaling)
    
    # Clamp peaks falling outside the range to the range
    fit_mu <- fit$par[1]
    if(fit_mu < 1){
      fit_mu <- 1
    } else if(fit_mu > conN){
      fit_mu <- conN
    }
    
    # Round values to discretise 
    fit_mu <- round(fit_mu)
    
    # Add back to the xii
    grad_xii$data$cortex_right[j]    <- fit_mu
    sigma_xii$data$cortex_right[j]   <- fit$par[2]
    scaling_xii$data$cortex_right[j] <- fit$par[3]
    value_xii$data$cortex_right[j]   <- fit$value
  }
  
  ## subcort
  for(j in 1:length(grad_xii$data$subcort)){
    # Initialise values
    values <- c()
    
    # Loop through all levels
    for(i in 1:conN){
      values[i] <- xii_list[[i]]$data$subcort[j]
    }
    
    # Make all values positive
    values <- values - min(values)
    
    # Fit Gaussian function
    fit <- fit_gaussian_curve(x, values, mu, sig, scaling)
    
    # Clamp peaks falling outside the range to the range
    fit_mu <- fit$par[1]
    if(fit_mu < 1){
      fit_mu <- 1
    } else if(fit_mu > conN){
      fit_mu <- conN
    }
    
    # Round values to discretise 
    fit_mu <- round(fit_mu)
    
    # Add back to the xii
    grad_xii$data$subcort[j]    <- fit_mu
    sigma_xii$data$subcort[j]   <- fit$par[2]
    scaling_xii$data$subcort[j] <- fit$par[3]
    value_xii$data$subcort[j]   <- fit$value
  }
  
  # Get the old labels
  old_key_colours <- grad_xii$meta$cifti$labels$`vertex areas`
  
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
  grad_xii$meta$cifti$labels$`vertex areas` <- new_key_colours
  
  # Save in list
  return(list(grad_xii = grad_xii, 
              sigma_xii = sigma_xii, 
              scaling_xii = scaling_xii,
              value_xii = value_xii))
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

bootstrap_analysis <- function(data, lm_formula, ddply_factors, nIter, colName, imageName){
  # Initialise results list
  results <- list()
  
  # Select correct column for analysis and create new data frame
  data$val     <- data[, colName]
  data2shuffle <- data
  
  # Calculate empirical values and save to list
  results$lm <- lm(lm_formula, data = data)
  numCoef    <- length(results$lm$coefficients) - 1 # Ignoring the intercept
  
  # Start cluster
  my.cluster <- parallel::makeCluster(detectCores() - 2, type = "PSOCK")
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  # Run parallel loop
  bootstrapped_values <- foreach(i = 1:nIter, .combine = 'c', .packages = 'plyr') %dopar% {
    data2shuffle$val <- sample(data$val)
    data2shuffle_agg <- ddply(data2shuffle, ddply_factors, summarise, val = mean(val))
    
    # Fit model
    temp_lm <- lm(lm_formula, data = data2shuffle_agg)
    
    # Add values 
    temp_est <- as.data.frame(matrix(as.numeric(temp_lm$coefficients)[-1], ncol = numCoef))
    names(temp_est) <- names(results$lm$coefficients)[-1]
    list(temp_est)
  }
  
  # Stop cluster again
  parallel::stopCluster(cl = my.cluster) 
  
  # Add to results
  results$bootstrapped_values <- as.data.frame(rbindlist(bootstrapped_values))
  
  # Save to disk  
  save(results, file = paste0("intermediate_data/", imageName))
  
  # Return value
  return(results)
}

# /* 
# ----------------------------- Calculate the gradients ---------------------------
# */
# Input for gradient analysis
conN       <- 6
folderName <- "SpaNov/OLMe_7T_SpaNov_gradient_6lvl_smo2_MSMAll"
copeNums   <- 8:13 # Because the copes are actually reversed relative to the 4 level version
fileName   <- "/stats/vwc/results_lvl2cope1_dat_ztstat_c1.dscalar.nii"
novLabels  <- paste0('lvl', 1:conN)
novFam_gradient      <- viridis(n = 6, option = "H", direction = -1)

# Load template file
# Folder where the images are
ciftiFolder <- "/SpaNov/OLMe_7T_SpaNov_gradient_6lvl_smo4_MSMAll/cope7.feat/stats/vwc/"

# Choose the files
zMap_file        <- paste0(path2imaging_results2, ciftiFolder, "results_lvl2cope1_dat_ztstat_c1.dscalar.nii")

# Load all maps as xiis
GLM1_zMap_xii <- read_cifti(zMap_file, brainstructures = "all", 
                            surfL_fname = surfLeft, 
                            surfR_fname = surfRight)

# Calculate the gradient map
against_avg  <- create_gradient_xiis(conN, folderName, copeNums, fileName, novLabels, "zstat", novFam_gradient)
against_avg2 <- create_gaussian_gradient_xiis(conN, folderName, copeNums, fileName, novLabels, "zstat", novFam_gradient)


# /* 
# ----------------------------- Write CIFTI files (wholebrain) ---------------------------
# */
# Write to disk as cifti
write_cifti(against_avg$grad_min_xii, cifti_fname = "cifti_results/SpaNov_gradient_wholebrain_min.dlabel.nii",
            surfL_fname = surfLeft,
            surfR_fname = surfRight,
            verbose = FALSE)
write_cifti(against_avg$grad_max_xii, cifti_fname = "cifti_results/SpaNov_gradient_wholebrain_max.dlabel.nii",
            surfL_fname = surfLeft,
            surfR_fname = surfRight,
            verbose = FALSE)

# Write to disk as .RData
save(against_avg, file = "cifti_results/SpaNovGradient_against_avg.RData")
save(against_avg2, file = "cifti_results/SpaNovGradient_against_avg_Gaussian.RData")

# /* 
# ----------------------------- Write CIFTI files (PMC) ---------------------------
# */
# Minimum
## Load the file
load("intermediate_data/SpaNov_cortical_gradient_analysis_min.RData")

## Write the L & R to disk
cifti_fname <- "cifti_results/SpaNovGradient_min_L_1.dlabel.nii"
write_cifti(xifti = gradient_xiftis_L[[1]], cifti_fname = cifti_fname)

cifti_fname <- "cifti_results/SpaNovGradient_min_R_1.dlabel.nii"
write_cifti(xifti = gradient_xiftis_R[[1]], cifti_fname = cifti_fname)

# Code to investigate other gradient
# for(i in 1:67){
#   # Print i
#   print(i)
#   view_cifti_surface(gradient_xiftis_R[[i]])
#   
#   # Wait for user input
#   readline(prompt="Press [enter] to continue")
# }

# Combine both to one CIFTI
new_grad_xii <- gradient_xiftis_R[[1]]
new_grad_xii$data$cortex_left <-  gradient_xiftis_L[[1]]$data$cortex_left

# We also want to add all vertices that are zero in Grad 1 but not 
# zero in Grad 25 to be replaced
bool_index <- new_grad_xii$data$cortex_right == 0 & gradient_xiftis_R[[34]]$data$cortex_right != 0
new_grad_xii$data$cortex_right[bool_index] <- gradient_xiftis_R[[34]]$data$cortex_right[bool_index]

# Now we will cut-off all unnecessary regions that are clearly not part of the 
# gradients. Especially the ones on the lateral surface of the brain
regions2exclude <- c("VIP", "LIPv", "MIP", "LIPd", "IP0", "IPS1", "V3B", "V3CD", 
                     "PGp", "PGi", "TPOJ3", "MT", "V4t", "LO2", "LO3", "V4", "V3A",
                     "LO1", "PFm", 'IP1', "IP2")

regions2include <- c("7Pm", "POS2", "POS1", "RSC", "7m", "PCV", "5mv", "23c", 
                     "23d", "31a", "d23ab", "v23ab", "31pv", "31pd", "DVT")

# I checked if all the names are written correctly using: 
# all( regions2exclude %in% parcel_names$V2)

# MMP IDs to remove and create boolean index
MMP_IDs    <- parcel_names$V1[parcel_names$V2 %in% regions2include]

# Left
bool_index <- MMP_xii$data$cortex_left %in% (MMP_IDs + 180)
new_grad_xii$data$cortex_left[!bool_index] <- 0

# Right
bool_index <- MMP_xii$data$cortex_right %in% (MMP_IDs)
new_grad_xii$data$cortex_right[!bool_index] <- 0

# Write new file
cifti_fname <- "cifti_results/SpaNov_gradient_PMC_min.dlabel.nii"
write_cifti(xifti = new_grad_xii, cifti_fname = cifti_fname)

# /* 
# ----------------------------- Extract & save hippocampal gradient ---------------------------
# */
# Get the hippocampal data
HC_data <- get_HC_data_from_xii(against_avg)

# Save
save(HC_data, file = "intermediate_data/SpaNov_gradient_data.RData")
