# This script runs the cortical gradient analysis by 
# 1. Finding clusters in a gradient map
# 2. Then by going from one vertex of each cluster to see through how many other cluster one travel
# Background links:
# https://www.humanconnectome.org/software/workbench-command/-surface-geodesic-distance
# https://www.nitrc.org/forum/message.php?msg_id=38001
# If a stack error are encountered then this might solve the issue:
# https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit
# Run this before script: ulimit -s 96000 # enlarge stack limit to 16 megs
# I only tested this script on my ubuntu laptop. For this to work, I had to change
# the hardlimit of the maximum stack limit

options(expressions = 500000)

# /* 
# ----------------------------- Parameters ---------------------------
# */
# Minimum cluster size to serve as a seed cluster
clusterCutOff <- 10

# /* 
# ----------------------------- Setting up ---------------------------
# */
# Setting seed
set.seed(20231208)

# Libs
library(viridis)
library(assortedRFunctions)
library(plyr)

# Load the ciftiTools package and point to the Connectome Workbench 
library(ciftiTools)

# Use correct location based on computers
if(Sys.info()[4] == "DESKTOP-335I26I"){
  # Work laptop
  ciftiTools.setOption("wb_path", "C:/Program Files/workbench-windows64-v1.5.0/workbench/bin_windows64")
  path2imaging_results <- "D:/Seafile/imaging_results" # Work laptop
} else if(Sys.info()[4] == 'DESKTOP-91CQCSQ') {
  ciftiTools.setOption("wb_path", "D:/workbench/bin_windows64")
  path2imaging_results <- "D:/imaging_results" # Work desktop computer
} else if(Sys.info()[4] == 'alex-Zenbook-UX3404VA-UX3404VA') {
  ciftiTools.setOption("wb_path", "/usr/bin/wb_command")
  path2imaging_results <- "/media/alex/shared/Seafile/imaging_results" # Work laptop but ubuntu
} else if(Sys.info()[4] == "greengoblin"){
  # Main desktop PC (Linux)
  ciftiTools.setOption("wb_path", "/usr/bin/wb_command") 
  path2imaging_results <- "/media/alex/work/Seafile/imaging_results" 
} else if(Sys.info()[4] == "GREEN-GOBLIN-WI"){
  # Main desktop PC (Linux)
  ciftiTools.setOption("wb_path", "C:/Program Files/workbench/bin_windows64")
  path2imaging_results <- "E:/Seafile/imaging_results" 
} else {
  # Personal laptop (Windows)
  ## Setting paths to workbench installation
  ciftiTools.setOption("wb_path", "D:/Program Files/workbench/bin_windows64")
  
  ## Path to the imaging data
  path2imaging_results <- "D:/OLM/imaging_results"
}

# Loading the cifti files etc.
surfLeft  <- "data/ignore_fMRI_version1/sourceFiles/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii"
surfRight <- "data/ignore_fMRI_version1/sourceFiles/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii"
load("cifti_results/SpaNovGradient_against_avg.RData")

# Choose which map
grad_xifti <-  against_avg$grad_max_xii

# /* 
# ----------------------------- Functions ---------------------------
# */
# Find all neighbours in a vertex
find_neighbours <- function(current_vert, faces){
  # To which faces that the current vertex belong?
  faces_bool <- faces == current_vert
  
  # Convert matrix to an index
  face_index <- rowMeans(faces_bool) > 0 
  
  # Get all vertices that names by finding the unique values
  neighbours <- unique(c(faces[face_index, ]))
  
  # Remove current vertex from it
  neighbours <- neighbours[neighbours != current_vert]
  
  # Return
  return(neighbours)
}

# Find all neighboring vertices with the same value value aka clusters
find_whole_cluster <- function(current_vert, verbose = TRUE){
  if(verbose){
    cat("\rVisited verts:", sum(visited_verts), "| Visiting now:", current_vert)
  } 
  
  # Return if this vertex is already visited
  if(visited_verts[current_vert] != 0){
    return()
  }
  
  # Find all neighbours of this vertex
  neighbours    <- find_neighbours(current_vert, faces)
  
  # Exclude all neighbours that don't have the same value. This automatically excludes medial wall
  neighbours <- neighbours[grad_vert_full[neighbours] == cluster_value]
  
  # Sort neighbours
  neighbours <- sort(neighbours)
  
  # Get the number of neighbours and return if zero
  numNeighbours <- length(neighbours)
  if(numNeighbours == 0){
    return()
  }
  
  # Also return if all neighbours have been assigned already
  if(all(!is.na(cluster[neighbours]))){
    return()
  }
  
  # Update visited vertices
  visited_verts[current_vert] <- visited_verts[current_vert] + 1
  visited_verts <<- visited_verts
  
  # Loop through all vertices
  for(i in 1:numNeighbours){
    current_vert <- neighbours[i]
    find_whole_cluster(current_vert, verbose)
  }
}

# /* 
# ----------------------------- Find all clusters (left) ---------------------------
# */
# Create empty cluster
cluster_xifti <- grad_xifti
cluster_xifti$data$cortex_left  <- matrix(0, nrow = 29696, ncol = 1)
cluster_xifti$data$cortex_right <- matrix(0, nrow = 29716, ncol = 1)
cluster_xifti$data$subcort      <- matrix(0, nrow = 31870, ncol = 1)

# Left hemisphere
verts <- cluster_xifti$surf$cortex_left$vertices
faces <- cluster_xifti$surf$cortex_left$faces

# Get the information from the gradient map
grad_vert    <- grad_xifti$data$cortex_left
medial_wall  <- grad_xifti$meta$cortex$medial_wall_mask$left
cluster      <- matrix(NA, nrow = nrow(verts), ncol = 1)
cluster_size <-  matrix(NA, nrow = nrow(verts), ncol = 1)

# Create matrix that includes the medial wall vertices
grad_vert_full <- matrix(NA, nrow = nrow(verts), ncol = 1)
grad_vert_full[medial_wall] <- grad_vert

# Start with the first vertex
cluster_vert    <- 1
cluster_counter <- 1

# Loop until no vertices is unassigned
while(any(is.na(cluster))){
  # Only try to run if the current vert is not part of the medial wall
  if(!is.na(grad_vert_full[cluster_vert])){
    # Create a new visited cluster matrix
    visited_verts <- matrix(0, nrow = nrow(verts), ncol = 1)
    cluster_value <- grad_vert_full[cluster_vert]

    # Find the cluster
    find_whole_cluster(cluster_vert)
    
    # Update cluster matrix and counter
    cat("\n\nCurr. number of clusters:", cluster_counter, 
        "\nCurr. vertex:", cluster_vert, 
        "\nSize of this cluster: ", sum(visited_verts != 0), 
        "\nVertices remaining:", sum(is.na(cluster)))
    cat("\n")
    
    # Check if any vertices were even visted
    if(sum(visited_verts != 0) > 0){
      cluster[visited_verts != 0] <- cluster_counter
      cluster_size[visited_verts != 0] <- sum(visited_verts != 0)
    } else {
      cluster[cluster_vert] <- cluster_counter
      cluster_size[cluster_vert] <- 1
    }
    cluster_counter <- cluster_counter + 1
    
  } else {
    cluster[cluster_vert] <- 0
    cluster_size[cluster_vert] <- 1
  }
  
  # Move to the next vertex that is not already used
  cluster_vert <- which(is.na(cluster))[1]
}


# /* 
# ----------------------------- Add results back to xifti (left) ---------------------------
# */
# Add data back
cluster_xifti$data$cortex_left <- matrix(cluster[medial_wall], ncol = 1)

# Get the maximum cluster number in left as this has to be added to the right hemisphere
# to get unique key values
max_cluster_left <- max(cluster_xifti$data$cortex_left)

# Change the colours
conN       <- length_uniq(cluster[medial_wall])
conColours <- sample(viridis(n = conN, option = "H"))
keyLabels  <- unique(cluster[medial_wall])

# Get the old labels
old_key_colours <- cluster_xifti$meta$cifti$labels$`vertex areas`

# Create new colours
new_key_colours <- old_key_colours[-(2:nrow(old_key_colours)),]
row.names(new_key_colours)[1] <- ""

for(i in 1:conN){
  # Set the colours
  RGB_Col <- col2rgb(conColours[i])/255
  temp_key_colour <- data.frame(Key = keyLabels[i],
                                Red = RGB_Col[1],
                                Green = RGB_Col[2],
                                Blue = RGB_Col[3],
                                Alpha = 1)

  new_key_colours <- rbind(new_key_colours, temp_key_colour)
}

# Add the labels
row.names(new_key_colours) <- c("???", paste0("L_", keyLabels))

# Add back to xifti
cluster_xifti$meta$cifti$labels$`vertex areas` <- new_key_colours

# /* 
# ----------------------------- Find all clusters (right) ---------------------------
# */
# Create empty cluster
# Right hemisphere
verts <- cluster_xifti$surf$cortex_right$vertices
faces <- cluster_xifti$surf$cortex_right$faces

# Get the information from the gradient map
grad_vert    <- grad_xifti$data$cortex_right
medial_wall  <- grad_xifti$meta$cortex$medial_wall_mask$right
cluster      <- matrix(NA, nrow = nrow(verts), ncol = 1)
cluster_size <-  matrix(NA, nrow = nrow(verts), ncol = 1)

# Create matrix that includes the medial wall vertices
grad_vert_full <- matrix(NA, nrow = nrow(verts), ncol = 1)
grad_vert_full[medial_wall] <- grad_vert

# Start with the first vertex
cluster_vert    <- 1
cluster_counter <- 1

# Loop until no vertices is unassigned
while(any(is.na(cluster))){
  # Only try to run if the current vert is not part of the medial wall
  if(!is.na(grad_vert_full[cluster_vert])){
    # Create a new visited cluster matrix
    visited_verts <- matrix(0, nrow = nrow(verts), ncol = 1)
    cluster_value <- grad_vert_full[cluster_vert]
    
    # Find the cluster
    find_whole_cluster(cluster_vert)
    
    # Update cluster matrix and counter
    cat("\n\nCurr. number of clusters:", cluster_counter, 
        "\nCurr. vertex:", cluster_vert, 
        "\nSize of this cluster: ", sum(visited_verts != 0), 
        "\nVertices remaining:", sum(is.na(cluster)))
    cat("\n")
    
    # Check if any vertices were even visted
    if(sum(visited_verts != 0) > 0){
      cluster[visited_verts != 0] <- cluster_counter
      cluster_size[visited_verts != 0] <- sum(visited_verts != 0)
    } else {
      cluster[cluster_vert] <- cluster_counter
      cluster_size[cluster_vert] <- 1
    }
    cluster_counter <- cluster_counter + 1
    
  } else {
    cluster[cluster_vert] <- 0
    cluster_size[cluster_vert] <- 1
  }
  
  # Move to the next vertex that is not already used
  cluster_vert <- which(is.na(cluster))[1]
}

# /*
# ----------------------------- Add results back to xifti (right) ---------------------------
# */
# Re-number so the key values are unique from the left hemisphere
cluster[cluster != 0] <- cluster[cluster != 0] + max_cluster_left

# Add data back
cluster_xifti$data$cortex_right <- matrix(cluster[medial_wall], ncol = 1)

# Change the colours
conN <- length_uniq(cluster[medial_wall])
conColours <- sample(viridis(n = conN, option = "H"))
keyLabels <- unique(cluster[medial_wall])

# Get the old labels
old_key_colours <- cluster_xifti$meta$cifti$labels$`vertex areas`

# Create new colours
new_key_colours <- old_key_colours[-(1:nrow(old_key_colours)),]

for(i in 1:conN){
  # Set the colours
  RGB_Col <- col2rgb(conColours[i])/255
  temp_key_colour <- data.frame(Key = keyLabels[i],
                                Red = RGB_Col[1],
                                Green = RGB_Col[2],
                                Blue = RGB_Col[3],
                                Alpha = 1)

  new_key_colours <- rbind(new_key_colours, temp_key_colour)
}

# Add the labels
row.names(new_key_colours) <- c(paste0("R_", keyLabels))

# Add to old labels that have already been set for the left hemisphere
new_key_colours <- rbind(old_key_colours, new_key_colours)

# Add back to xifti
cluster_xifti$meta$cifti$labels$`vertex areas` <- new_key_colours

# /*
# ----------------------------- Save the temp. results ---------------------------
# */
# Write CIFTI
write_cifti(xifti = cluster_xifti, 
            cifti_fname = "cifti_results/SpaNovGradient_against_avg_max_clusters.dlabel.nii")

# /*
# ----------------------------- Find gradient function ---------------------------
# */
# This function is nearly identical to the find_whole_cluster function but with 
# slight changes that allow going to neighbours that have different lvls as long
# as the lvls are N - 1. One exception is for Lvl 5 which is skipped in the PMC.
find_gradient <- function(current_vert, verbose = TRUE){
  if(verbose){
    cat("\rVisited verts:", sum(visited_verts), "| Visiting now:", current_vert)
  } 
  
  # Return if this vertex is already visited
  if(visited_verts[current_vert] != 0){
    return()
  }
  
  # Find all neighbours of this vertex
  neighbours    <- find_neighbours(current_vert, faces)
  
  # Get the current lvl of this vertex
  current_lvl <- grad_vert_full[current_vert]
  target_lvl  <- current_lvl - 1
  # 5 is skipped for the main gradient, we therefore allow transitions from 6 to 5
  # as well 6 directly
  if(target_lvl == 5){
    target_lvl <- c(4, 5)
  } else if(target_lvl < 1){
    target_lvl <- 1
  }
  permissible_lvls <- c(current_lvl, target_lvl)
  
  # Get lvls of the neighbours 
  neighbours_lvl <- grad_vert_full[neighbours]
  
  # Subset to permissible neighbours
  neighbours <- neighbours[neighbours_lvl %in% permissible_lvls]
  
  # Sort neighbours
  neighbours <- sort(neighbours)
  
  # Get the number of neighbours and return if zero
  numNeighbours <- length(neighbours)
  if(numNeighbours == 0){
    return()
  }
  
  # Update visited vertices
  visited_verts[current_vert] <- visited_verts[current_vert] + 1
  visited_verts <<- visited_verts
  
  # Loop through all vertices
  for(i in 1:numNeighbours){
    current_vert <- neighbours[i]
    find_gradient(current_vert, verbose)
  }
}

# /*
# ----------------------------- Running gradient analysis for each lvl 6 cluster left ---------------------------
# */
# Get one vertex for each Level 6 cluster
## Get all information 
verts        <- cluster_xifti$surf$cortex_left$vertices
faces        <- cluster_xifti$surf$cortex_left$faces
cluster_vert <- cluster_xifti$data$cortex_left 
grad_vert    <- grad_xifti$data$cortex_left
medial_wall  <- grad_xifti$meta$cortex$medial_wall_mask$left
cluster_vert_full              <- matrix(NA, nrow = nrow(verts), ncol = 1)
cluster_vert_full[medial_wall] <- cluster_vert
grad_vert_full                 <- matrix(NA, nrow = nrow(verts), ncol = 1)
grad_vert_full[medial_wall]    <- grad_vert

## Create a df
cluster_df_L      <- data.frame(cluster = cluster_vert_full, lvl = grad_vert_full)
cluster_df_L$vert <- 1:nrow(cluster_df_L) 

## Remove everything unnecessary and get the vertex per cluster
cluster_df_L  <- cluster_df_L[!is.na(cluster_df_L$cluster) & cluster_df_L$lvl == 6, ]
cluster_df_L  <- ddply(cluster_df_L, c("cluster", "lvl"), summarise, vert = vert[1], size = length(lvl))

# Remove clusters that are too small 
cluster_df_L    <- cluster_df_L[cluster_df_L$size >= clusterCutOff, ]

# List for the xifti vars and prepate cluster_df_L for extra info
gradient_xiftis_L          <- list()
cluster_df_L$num_lvls      <- NA
cluster_df_L$gradient_size <- NA # Use two methods to calculate size to check algorithm
cluster_df_L$visited_sum   <- NA # Use two methods to calculate size to check algorithm

# Print message
cat("\nLooking for gradients in the left hemisphere...\n")

# Start time for progress bar
startTime <- Sys.time()

# Loop through all remaining clusters
for(i in 1:nrow(cluster_df_L)){
  # Reset the visited vert matrix
  visited_verts <- matrix(0, nrow = nrow(verts), ncol = 1)
  
  # Run the algorithm
  success <- tryCatch({
    find_gradient(cluster_df_L$vert[i], verbose = FALSE)
    TRUE
  }, error = function(error_condition) {
    cat(as.character(error_condition))
    FALSE
  })
  
  # Process the results from the algorithm
  if(success){
    ## Get all levels that are part of this gradient
    gradient_lvls       <- grad_vert_full[visited_verts != 0]
    num_gradient_levels <- length_uniq(gradient_lvls)
    visited_sum         <- sum(visited_verts) # Use two methods to calculate size to check algorithm
    gradient_size       <- sum(visited_verts != 0) # Use two methods to calculate size to check algorithm
    
    ## Add information back to cluster_df_L
    cluster_df_L$num_lvls[i]      <- num_gradient_levels
    cluster_df_L$gradient_size[i] <- gradient_size
    cluster_df_L$visited_sum[i]   <- visited_sum
    
    ## Create a new xifti var for this gradient
    gradient_xiftis_L[[i]] <- grad_xifti
    gradient_xiftis_L[[i]]$data$cortex_left[visited_verts[medial_wall] == 0] <- 0
    gradient_xiftis_L[[i]]$data$cortex_right <- matrix(0, nrow = 29716, ncol = 1)
  } else {
    cluster_df_L$num_lvls[i]      <- 0
    cluster_df_L$gradient_size[i] <- NA
    cluster_df_L$visited_sum[i]   <- NA
    
    ## Create a new xifti var for this gradient
    gradient_xiftis_L[[i]] <- NA
  }
  
  # Print progress
  progressDisplay(i, nrow(cluster_df_L), startTime)
}

# New line
cat("\n")


# /*
# ----------------------------- Colour the seed clusters ---------------------------
# */
seed_clusters <- cluster_xifti
seed_clusters$data$cortex_left <- matrix(0, nrow = 29696, ncol = 1)
seed_clusters$data$cortex_right <- matrix(0, nrow = 29716, ncol = 1)
seed_clusters$data$subcort <- matrix(0, nrow = 31870, ncol = 1)

# Loop through all clusters
for(i in 1:nrow(cluster_df_L)){
  cluster_index <- cluster_xifti$data$cortex_left == cluster_df_L$cluster[i]
  seed_clusters$data$cortex_left[cluster_index] <- cluster_df_L$num_lvls[i]
}

# Change the colours
conN <- 6
conColours <- sample(viridis(n = conN, option = "H"))
keyLabels <- 1:6

# Get the old labels
old_key_colours <- cluster_xifti$meta$cifti$labels$`vertex areas`

# Create new colours
new_key_colours <- old_key_colours[-(2:nrow(old_key_colours)),]

for(i in 1:conN){
  # Set the colours
  RGB_Col <- col2rgb(conColours[i])/255
  temp_key_colour <- data.frame(Key = keyLabels[i],
                                Red = RGB_Col[1],
                                Green = RGB_Col[2],
                                Blue = RGB_Col[3],
                                Alpha = 1)
  
  new_key_colours <- rbind(new_key_colours, temp_key_colour)
}

row.names(new_key_colours) <- c("???", paste0("Num. lvls = ", keyLabels))

# Add back to xifti
seed_clusters$meta$cifti$labels$`vertex areas` <- new_key_colours

# /*
# ----------------------------- Running gradient analysis for each lvl 6 cluster right ---------------------------
# */
# Get one vertex for each Level 6 cluster
## Get all information 
verts        <- cluster_xifti$surf$cortex_right$vertices
faces        <- cluster_xifti$surf$cortex_right$faces
cluster_vert <- cluster_xifti$data$cortex_right 
grad_vert    <- grad_xifti$data$cortex_right
medial_wall  <- grad_xifti$meta$cortex$medial_wall_mask$right
cluster_vert_full              <- matrix(NA, nrow = nrow(verts), ncol = 1)
cluster_vert_full[medial_wall] <- cluster_vert
grad_vert_full                 <- matrix(NA, nrow = nrow(verts), ncol = 1)
grad_vert_full[medial_wall]    <- grad_vert

## Create a df
cluster_df_R      <- data.frame(cluster = cluster_vert_full, lvl = grad_vert_full)
cluster_df_R$vert <- 1:nrow(cluster_df_R) 

## Remove everything unnecessary and get the vertex per cluster
cluster_df_R  <- cluster_df_R[!is.na(cluster_df_R$cluster) & cluster_df_R$lvl == 6, ]
cluster_df_R  <- ddply(cluster_df_R, c("cluster", "lvl"), summarise, vert = vert[1], size = length(lvl))

# Cut-off for minimum cluster size
cluster_df_R    <- cluster_df_R[cluster_df_R$size >= clusterCutOff, ]

# List for the xifti vars and prepate cluster_df_R for extra info
gradient_xiftis_R          <- list()
cluster_df_R$num_lvls      <- NA
cluster_df_R$gradient_size <- NA # Use two methods to calculate size to check algorithm
cluster_df_R$visited_sum   <- NA # Use two methods to calculate size to check algorithm

# Print message
cat("\nLooking for gradients in the right hemisphere...\n")

# Start time for progress bar
startTime <- Sys.time()

# Loop through all remaining clusters
for(i in 1:nrow(cluster_df_R)){
  # Reset the visited vert matrix
  visited_verts <- matrix(0, nrow = nrow(verts), ncol = 1)
  
  # Run the algorithm
  success <- tryCatch({
    find_gradient(cluster_df_R$vert[i], verbose = FALSE)
    TRUE
  }, error = function(error_condition) {
    cat(as.character(error_condition))
    FALSE
  })
  
  # Process the results from the algorithm
  if(success){
    ## Get all levels that are part of this gradient
    gradient_lvls       <- grad_vert_full[visited_verts != 0]
    num_gradient_levels <- length_uniq(gradient_lvls)
    visited_sum         <- sum(visited_verts) # Use two methods to calculate size to check algorithm
    gradient_size       <- sum(visited_verts != 0) # Use two methods to calculate size to check algorithm
    
    ## Add information back to cluster_df_L
    cluster_df_R$num_lvls[i]      <- num_gradient_levels
    cluster_df_R$gradient_size[i] <- gradient_size
    cluster_df_R$visited_sum[i]   <- visited_sum
    
    ## Create a new xifti var for this gradient
    gradient_xiftis_R[[i]] <- grad_xifti
    gradient_xiftis_R[[i]]$data$cortex_left <- matrix(0, nrow = 29696, ncol = 1)
    gradient_xiftis_R[[i]]$data$cortex_right[visited_verts[medial_wall] == 0] <- 0
  } else {
    cluster_df_R$num_lvls[i]      <- 0
    cluster_df_R$gradient_size[i] <- NA
    cluster_df_R$visited_sum[i]   <- NA
    
    ## Create a new xifti var for this gradient
    gradient_xiftis_R[[i]] <- NA
  }
  
  # Print progress
  progressDisplay(i, nrow(cluster_df_R), startTime)
}

# New line
cat("\n")


# /*
# ----------------------------- Colour the seed clusters ---------------------------
# */
# Loop through all clusters
for(i in 1:nrow(cluster_df_R)){
  cluster_index <- cluster_xifti$data$cortex_right == cluster_df_R$cluster[i]
  seed_clusters$data$cortex_right[cluster_index] <- cluster_df_R$num_lvls[i]
}

# /*
# ----------------------------- Save results ---------------------------
# */
# Create .RData image
save.image("intermediate_data/SpaNov_cortical_gradient_analysis_max.RData")

write_cifti(xifti = seed_clusters, 
            cifti_fname = "cifti_results/SpaNovGradient_against_avg_max_seed_clusters.dlabel.nii")