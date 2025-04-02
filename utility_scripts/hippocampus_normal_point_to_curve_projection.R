# This script projects each voxel on to middle plane/line of the hippocampus. 

# Notes:
# For what ever reason the th right hippocampus point 84 gives are weird value that
# is not 90 degrees from the tangent even if there is such value. The easiest way
# to fix this is issue that I came up with is to just change the starting value 
# just for this point only.

# /* 
# ----------------------------- Libraries and parameter ---------------------------
# */
# Set seed
set.seed(19911225)

# Fitting parameters
stepSize  <- 0.0001

# Colours 
curveColour      <- "black"
unfinishedColour <- "#fc2803"
beingFitColour   <- "#ffff66"
finishedColour   <- "#99ff66"
porjectedPoints  <- "grey"
lineSize  <- 1
pointSize <- 2

# Libraries
library(ggplot2)
library(assortedRFunctions)
library(cowplot)
library(stringr)
library(plyr)
library(mgcv)

# /* 
# ----------------------------- Load hippocampus data---------------------------
# */
# Load the extracted MNI coordinates
MNI_coord <- read.csv("otherStuff/cifti_subcortical_MNI152_coordinates.csv")

# Get hippocampal values
HC_L <- MNI_coord[MNI_coord$region == "Hippocampus-L", ]
HC_R <- MNI_coord[MNI_coord$region == "Hippocampus-R", ]

# /* 
# ----------------------------- Functions for the algorithm and others ---------------------------
# */
plot_function <- function(){
  # Scatter plot
  p1 <- ggplot() + 
    geom_line(data = sampleCurve, mapping = aes(x = y, y = z_prime), colour = curveColour, size = lineSize) +
    geom_point(data = todo_df, mapping = aes(x = y, y = z), fill = unfinishedColour, pch = 21, size = pointSize) +
    #geom_text(data = unique_points, mapping = aes(x = y, y = z, label = id)) + 
    geom_point(data = finished_df, mapping = aes(x = y, y = z), fill = finishedColour, pch = 21, size = pointSize) +
    geom_point(data = projected_DF , mapping = aes(x = y, y = z_prime), fill = porjectedPoints, pch = 21, size = pointSize) +
    geom_segment(data = data.frame(), mapping = aes(x = currentPoint$y, 
                                                    y = currentPoint$z, 
                                                    xend = currentValue_vec[1], 
                                                    yend = currentValue_vec[2]), size = lineSize) +
    geom_point(data = currentPoint, mapping = aes(x = y, y = z), fill = beingFitColour, pch = 21, size = pointSize) +
    geom_abline(slope = mb[1], intercept = mb[2], size = lineSize) +
    theme_void() +
    labs(title = "", x = "", y = "") +
    coord_fixed()
  
  # Circle plot
  p2 <- ggplot() + 
    geom_path(data = circleDat, mapping = aes(x = x, y = y), size = lineSize) +
    geom_segment(data = data.frame(), mapping = aes(x = origin_vec[1], 
                                                    y = origin_vec[2], 
                                                    xend = data_vec[1], 
                                                    yend = data_vec[2]), size = lineSize) +
    geom_segment(data = data.frame(), mapping = aes(x = origin_vec[1], 
                                                    y = origin_vec[2], 
                                                    xend = tan_vec2[1], 
                                                    yend = tan_vec2[2]), size = lineSize) +
    geom_point(data = vector_df, mapping = aes(x = x, y = y, colour = point), size = pointSize * 1.5) +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_void() +
    labs(title = "", x = "", y = "") +
    theme(legend.position = "none")
  
  print(plot_grid(p1, p2, ncol = 2))
}

calculate_slope_and_intercept <- function(x, y){
  # Calculate tangent & intercept
  m <- (y[2]-y[1])/(x[2]-x[1])
  b <- (m*x[1]-y[1])*-1
  
  return(c(m, b))
}

get_valuesOfTangent <- function(currentValue, stepSize, curve_model){
  y       <- currentValue[1] - stepSize
  z_prime <- predict(curve_model, newdata = data.frame(y = y))
  
  return(as.numeric(c(y, z_prime)))
}

# Draw a circle so the boundary can be plotted. 
circleDat <- circleFun(c(0,-0), 2, npoints = 100)

# Make sure that the circle is actually closing
circleDat <- rbind(circleDat, circleDat[2, ])

# Optimisation function
optim_fun <- function(data, par){
  # Get new value
  new_z_prime <- as.numeric(predict(curve_model, newdata = data.frame(y = par)))
  currentValue_vec    <- c(par, new_z_prime)
  
  tan_vec1       <- get_valuesOfTangent(currentValue_vec, stepSize, curve_model) # Vector for tangent
  data_vec       <- c(data$y, data$z) # Vector for data
  origin_vec     <- currentValue_vec # Same as startValue_vec but will be changed
  
  # Make all vectors relative to origin
  tan_vec2    <- tan_vec1 - origin_vec
  data_vec    <- data_vec - origin_vec
  origin_vec  <- origin_vec - origin_vec
  
  # Make the non-origin vectors unit length
  tan_vec2 <- unitVector(tan_vec2)
  data_vec <- unitVector(data_vec)
  
  # Calculate difference
  angle1             <- getAngleInDegreesFromPoint(tan_vec2[1], tan_vec2[2])
  angle2             <- getAngleInDegreesFromPoint(data_vec[1], data_vec[2])
  angleDiff          <- angularDifference(angle1, angle2)
  angleDiff_abs      <- abs(angleDiff)
  abs(angleDiff_abs - 90)^2
}

# /* 
# ----------------------------- Run projection for left hippocampus ---------------------------
# */
# Fit curve model as GAM to get the hippocampal plane following the overall shape
curve_model <- gam(formula = z ~ s(y, bs = "cs"), data = HC_L)

# Unique y and z points
unique_pairs  <- paste0(HC_L$y, "_", HC_L$z)
unique_points <- HC_L[!duplicated(unique_pairs), ]
unique_points$id <- 1:nrow(unique_points)

# Points to project
points2project <- sort(unique(unique_points$y))

# Predict y values 
minValue <- min(points2project)
maxValue <- max(points2project)
yLength  <- maxValue - minValue
sampleCurve         <- data.frame(y = seq(from = minValue, to = maxValue, length.out = 1000))
sampleCurve$z_prime <- predict(curve_model, newdata = sampleCurve)

# Data frames
projected_DF  <- data.frame(matrix(ncol = 2, nrow = 0))
names(projected_DF) <- c("y", "z_prime")
todo_df <- data.frame(y = unique_points$y, z = unique_points$z)
finished_df   <- data.frame(matrix(ncol = 2, nrow = 0))
names(finished_df) <- c("y", "z")

# Save angles and differences
endAngle     <- c()
end_y        <- c()
end_z  <- c()
endDiff      <- c()
converge     <- c()
converge_msg <- c()

# Loop through all points
for(i in 1:nrow(unique_points)){
  # Current point
  currentPoint <- data.frame(y = unique_points$y[i], z = unique_points$z[i], z_prime = NA)
  currentPoint$z_prime <- as.numeric(predict(curve_model, newdata = data.frame(y = currentPoint$y)))
  
  # Optimisation
  optimResult <- optim(par = currentPoint$y,
        fn = optim_fun,
        method = "L-BFGS-B",
        lower = minValue,
        upper = maxValue,
        data = currentPoint)
  
  # Update DFs for visualisations
  z_prime          <- as.numeric(predict(curve_model, newdata = data.frame(y = optimResult$par)))
  currentValue_vec <- c(optimResult$par, z_prime)
  tan_vec1         <- get_valuesOfTangent(currentValue_vec, stepSize, curve_model) # Vector for tangent
  data_vec         <- c(currentPoint$y, currentPoint$z) # Vector for data
  origin_vec       <- currentValue_vec # Same as startValue_vec but will be changed
  
  # Make all vectors relative to origin
  tan_vec2     <- tan_vec1 - origin_vec
  data_vec    <- data_vec - origin_vec
  origin_vec  <- origin_vec - origin_vec
  
  # Make the non-origion vectors unit lenght
  tan_vec2 <- unitVector(tan_vec2)
  data_vec <- unitVector(data_vec)
  
  # Create vector DF
  vector_df <- data.frame(x = c(origin_vec[1], tan_vec2[1], data_vec[1]),
                          y = c(origin_vec[2], tan_vec2[2], data_vec[2]),
                          point = c("origin", "tangent", "data"))
  
  # Calculate angles
  angle1             <- getAngleInDegreesFromPoint(tan_vec2[1], tan_vec2[2])
  angle2             <- getAngleInDegreesFromPoint(data_vec[1], data_vec[2])
  angleDiff          <- angularDifference(angle1, angle2)
  angleDiff_abs      <- abs(angleDiff)
  currentDiff        <- abs(angleDiff_abs - 90)
  currentDiff_signed <- angleDiff - 90

  # Calculate slope and intercept of tangent line
  mb <- calculate_slope_and_intercept(c(currentValue_vec[1], tan_vec1[1]),
                                      c(currentValue_vec[2], tan_vec1[2]))
  
  # Save results
  endAngle[i]     <- angleDiff_abs
  end_y[i]        <- optimResult$par
  end_z[i]        <- z_prime
  endDiff[i]      <- optimResult$value
  converge[i]     <- optimResult$convergence
  converge_msg[i] <- optimResult$message
  
  # Add to DF
  projected_DF <- rbind(projected_DF, data.frame(y = end_y[i], z_prime = end_z[i]))
  
  # Move row from unfinished to finish
  finished_df <- rbind(finished_df, todo_df[i, ])
  
  # Show current projection
  plot_function()
}


# Add back to points
unique_points$end_z <- end_z
unique_points$end_y <- end_y

# Sort by end_y
unique_points <- unique_points[order(unique_points$end_y, decreasing = TRUE),]

# Create new column/variable
unique_points$dist2previous <- 0

# Loop through the points
for(i in 2:nrow(unique_points)){
  previous <- unique_points[i - 1, c("end_y", "end_z")]
  current <- unique_points[i, c("end_y", "end_z")]
  
  # Calculate distance
  unique_points$dist2previous[i] <- dist(matrix(c(previous, current), ncol = 2, byrow = TRUE))
}

# Calculate position on plane by using cumsum
unique_points$position <-  cumsum(unique_points$dist2previous)

# Create plot
p1 <- ggplot(unique_points, aes(x = y, y = z)) + 
  geom_rect(aes(xmin =  y - 1, xmax = y + 1, ymin = z - 1, ymax = z + 1, fill = position)) +
  geom_line(data = sampleCurve, mapping = aes(x = y, y = z_prime), colour = "white", size = 1) +
  scale_fill_viridis_c(option = "C", direction = 1) +
  coord_fixed() +
  labs(title = "Projection along the \nshape (left)", x = "MNI y", y = "MNI z") +
  theme_void() +
  theme(legend.position = "none")

p2 <- ggplot(unique_points, aes(x = y, y = z)) + 
  geom_rect(aes(xmin =  y - 1, xmax = y + 1, ymin = z - 1, ymax = z + 1, fill = y)) +
  geom_line(data = sampleCurve, mapping = aes(x = y, y = z_prime), colour = "white", size = 1) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  coord_fixed() +
  labs(title = "MNI-y \ncoordinate (left)", x = "MNI y", y = "MNI z") +
  theme_void() +
  theme(legend.position = "none")

p3 <- ggplot(unique_points, aes(x = y, y = z)) + 
  geom_rect(aes(xmin =  y - 1, xmax = y + 1, ymin = z - 1, ymax = z + 1, fill = z)) +
  geom_line(data = sampleCurve, mapping = aes(x = y, y = z_prime), colour = "white", size = 1) +
  scale_fill_viridis_c(option = "C", direction = 1) +
  coord_fixed() +
  labs(title = "MNI-z \ncoordinate (left)", x = "MNI y", y = "MNI z") +
  theme_void() +
  theme(legend.position = "none")

# Save unique_points for this hemisphere
unique_points_L <- unique_points

# /* 
# ----------------------------- Run projection for right hippocampus ---------------------------
# */
# Fit curve model as GAM to get the hippocampal plane following the overall shape
curve_model <- gam(formula = z ~ s(y, bs = "cs"), data = HC_R)

# Unique y and z points
unique_pairs  <- paste0(HC_R$y, "_", HC_R$z)
unique_points <- HC_R[!duplicated(unique_pairs), ]
unique_points$id <- 1:nrow(unique_points)

# Points to project
points2project <- sort(unique(unique_points$y))

# Predict y values 
minValue <- min(points2project)
maxValue <- max(points2project)
yLength  <- maxValue - minValue
sampleCurve <- data.frame(y = seq(from = minValue, to = maxValue, length.out = 1000))
sampleCurve$z_prime <- predict(curve_model, newdata = sampleCurve)

# Data frames
projected_DF  <- data.frame(matrix(ncol = 2, nrow = 0))
names(projected_DF) <- c("y", "z_prime")
todo_df <- data.frame(y = unique_points$y, z = unique_points$z)
finished_df   <- data.frame(matrix(ncol = 2, nrow = 0))
names(finished_df) <- c("y", "z")

# Save angles and differences
endAngle     <- c()
end_y        <- c()
end_z        <- c()
endDiff      <- c()
converge     <- c()
converge_msg <- c()

# Loop through all points
for(i in 1:nrow(unique_points)){
  # Current point
  currentPoint <- data.frame(y = unique_points$y[i], z = unique_points$z[i], z_prime = NA)
  currentPoint$z_prime <- as.numeric(predict(curve_model, newdata = data.frame(y = currentPoint$y)))
  
  # Optimisation
  
  optimResult <- optim(par = ifelse(i == 84, currentPoint$y + 1, currentPoint$y),
                       fn = optim_fun,
                       method = "L-BFGS-B",
                       lower = minValue,
                       upper = maxValue,
                       data = currentPoint)
  
  # Update DFs for visualisations
  z_prime  <- as.numeric(predict(curve_model, newdata = data.frame(y = optimResult$par)))
  currentValue_vec <- c(optimResult$par, z_prime)
  tan_vec1         <- get_valuesOfTangent(currentValue_vec, stepSize, curve_model) # Vector for tangent
  data_vec         <- c(currentPoint$y, currentPoint$z) # Vector for data
  origin_vec       <- currentValue_vec # Same as startValue_vec but will be changed
  
  # Make all vectors relative to origin
  tan_vec2     <- tan_vec1 - origin_vec
  data_vec    <- data_vec - origin_vec
  origin_vec  <- origin_vec - origin_vec
  
  # Make the non-origion vectors unit lenght
  tan_vec2 <- unitVector(tan_vec2)
  data_vec <- unitVector(data_vec)
  
  # Create vector DF
  vector_df <- data.frame(x = c(origin_vec[1], tan_vec2[1], data_vec[1]),
                          y = c(origin_vec[2], tan_vec2[2], data_vec[2]),
                          point = c("origin", "tangent", "data"))
  
  # Calculate angles
  angle1             <- getAngleInDegreesFromPoint(tan_vec2[1], tan_vec2[2])
  angle2             <- getAngleInDegreesFromPoint(data_vec[1], data_vec[2])
  angleDiff          <- angularDifference(angle1, angle2)
  angleDiff_abs      <- abs(angleDiff)
  currentDiff        <- abs(angleDiff_abs - 90)
  currentDiff_signed <- angleDiff - 90
  
  # Calculate slope and intercept of tangent line
  mb <- calculate_slope_and_intercept(c(currentValue_vec[1], tan_vec1[1]),
                                      c(currentValue_vec[2], tan_vec1[2]))
  
  # Save results
  endAngle[i]     <- angleDiff_abs
  end_y[i]        <- optimResult$par
  end_z[i]        <- z_prime
  endDiff[i]      <- optimResult$value
  converge[i]     <- optimResult$convergence
  converge_msg[i] <- optimResult$message
  
  # Add to DF
  projected_DF <- rbind(projected_DF, data.frame(y = end_y[i], z_prime = end_z[i]))
  
  # Move row from unfinished to finish
  finished_df <- rbind(finished_df, todo_df[i, ])
  
  # Show current projection
  plot_function()
}

# Add back to points
unique_points$end_z <- end_z
unique_points$end_y <- end_y

# Sort by end_y
unique_points <- unique_points[order(unique_points$end_y, decreasing = TRUE),]

# Create new column/variable
unique_points$dist2previous <- 0

# Loop through the points
for(i in 2:nrow(unique_points)){
  previous <- unique_points[i - 1, c("end_y", "end_z")]
  current <- unique_points[i, c("end_y", "end_z")]
  
  # Calculate distance
  unique_points$dist2previous[i] <- dist(matrix(c(previous, current), ncol = 2, byrow = TRUE))
}

# Calculate position on plane by using cumsum
unique_points$position <-  cumsum(unique_points$dist2previous)

# Create plot
p4 <- ggplot(unique_points, aes(x = y, y = z)) + 
  geom_rect(aes(xmin =  y - 1, xmax = y + 1, ymin = z - 1, ymax = z + 1, fill = position)) +
  geom_line(data = sampleCurve, mapping = aes(x = y, y = z_prime), colour = "white", size = 1) +
  scale_fill_viridis_c(option = "C", direction = 1) +
  coord_fixed() +
  labs(title = "Projection along the \nshape (right)", x = "MNI y", y = "MNI z")+
  theme(legend.position = "none")

p5 <- ggplot(unique_points, aes(x = y, y = z)) + 
  geom_rect(aes(xmin =  y - 1, xmax = y + 1, ymin = z - 1, ymax = z + 1, fill = y)) +
  geom_line(data = sampleCurve, mapping = aes(x = y, y = z_prime), colour = "white", size = 1) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  coord_fixed() +
  labs(title = "MNI-y \ncoordinate (right)", x = "MNI y", y = "MNI z")+
  theme(legend.position = "none")

p6 <- ggplot(unique_points, aes(x = y, y = z)) + 
  geom_rect(aes(xmin =  y - 1, xmax = y + 1, ymin = z - 1, ymax = z + 1, fill = z)) +
  geom_line(data = sampleCurve, mapping = aes(x = y, y = z_prime), colour = "white", size = 1) +
  scale_fill_viridis_c(option = "C", direction = 1) +
  coord_fixed() +
  labs(title = "MNI-z \ncoordinate (right)", x = "MNI y", y = "MNI z")+
  theme(legend.position = "none")

# Save unique_points for this hemisphere
unique_points_R <- unique_points

# /* 
# ----------------------------- Visualise and save ---------------------------
# */
# Visualise the difference
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)


# Combine results and save
projected_HC <- rbind(unique_points_R, unique_points_L)
save(projected_HC, file = "otherStuff/projected_HC.RData")