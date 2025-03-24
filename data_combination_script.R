# This script loads all the data that is combined for the OLM project fMRI task (version 1) and saves this for sharing & analysis.
# The data is anonymised by removing identifiable information
# /* 
# ----------------------------- Preparations ---------------------------
# */
# VVVVVVVVVVVVVVVVVVVVVVVVV Needs specification VVVVVVVVVVVVVVVVVVVVVVVV
# Below are paths that need to be changed in order to run this script.
# Currently, it's only set-up to run on varius computers I use.

# Use correct location based on computer
if(Sys.info()[4] == "DESKTOP-335I26I"){
  # Work laptop - Windows
  setwd("D:/researchProjects/OLM_project/analysis")
  path2input_questions <- 'D:/researchProjects/OLM_project/preparation/OLM_versions/standard_fMRI/version1/session1/trialData/inputFiles/1_questionnaires/QP_input.txt'
} else if(Sys.info()[4] == 'DESKTOP-91CQCSQ') {
  # Work desktop
} else if(Sys.info()[4] == 'alex-Zenbook-UX3404VA-UX3404VA') {
  # Work laptop - Ubuntu
  # Specify analysis folder on the OLM_project as WD
  setwd("/media/alex/shared/researchProjects/OLM_project/analysis")
  path2input_questions <- '/media/alex/shared/researchProjects/OLM_project/preparation/OLM_versions/standard_fMRI/version1/session1/trialData/inputFiles/1_questionnaires/QP_input.txt'
} else {
  # Paths to the input .csv for the questionnaire
  #path2input_questions <- 'C:/Users/admin/Documents/GitHub/OLM_project/preparation/OLM_versions/standard_fMRI/version1/session1/trialData/inputFiles/1_questionnaires/QP_input.txt'
  #path2input_questions <- 'C:/Users/Alex/Documents/Work/researchProjects/OLM_project/preparation/OLM_versions/standard_fMRI/version1/session1/trialData/inputFiles/1_questionnaires/QP_input.txt'
}

# Paths that shouldn't need changing
## Path to where the combined data is saved
path2output <- "data/ignore_fMRI_version1/combined_data/" 

# Path where the data is saved
path2data <- "data/ignore_fMRI_version1/"

### ___________________________________________________________________
# Libs
# Install via devtools::install_github("JAQuent/assortedRFunctions", upgrade = "never")
library(assortedRFunctions)
library(jsonlite)
library(readxl)

# Load lookup table that contains the assignment of anonymised labels to the 
# participant information.
lookupTable  <- read.csv(paste0(path2data, "lookUpTable.csv"))

# Set locale
Sys.setlocale(category = "LC_ALL", locale = "chs")
###############

# /* 
# ----------------------------- Demographics ---------------------------
# */
# Load the spreadsheet that contains the demographic data
spreadsheet <- read_excel(paste0(path2data, "VP00157_PPT_Info_Final.xlsx"))
#spreadsheet <- spreadsheet[-1,] # Remove first row because this is only an example

# Subset with only the subjects that are used
spreadsheet <- spreadsheet[spreadsheet$`ID Number` %in% lookupTable$Rnum, ]

# Add anonymised subject code
spreadsheet$subject <- NA

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:nrow(lookupTable)){
  spreadsheet$subject[spreadsheet$`P Number` == lookupTable$HCnum[i]] <- lookupTable$anonKey[lookupTable$HCnum == lookupTable$HCnum[i]]
}

# Extract important information from the spreadsheet with regard to demographics
demographics <- data.frame(subject = spreadsheet$subject,
                           age = round(spreadsheet$Age, 1),
                           gender = spreadsheet$Gender,
                           timeBetween_sessions = spreadsheet$`7T MRI Scanning Date` - spreadsheet$`3T MRI Scanning Date`,
                           weekday_3T = weekdays(spreadsheet$`3T MRI Scanning Date`),
                           weekday_7T = weekdays(spreadsheet$`7T MRI Scanning Date`))

# Because the locale was set to Chinese the weekdays are also in Chinese which will be corrected
days_chinese_english <- data.frame(chinese = c("星期一", "星期二", "星期三", "星期四", "星期五", "星期六", "星期天"),
                                   english = c("monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday"))
# Loop through the days
for(i in 1:nrow(days_chinese_english)){
  # Translate the days
  demographics[demographics$weekday_3T == days_chinese_english$chinese[i] & !is.na(demographics$weekday_3T), 'weekday_3T'] <- days_chinese_english$english[i]
  demographics[demographics$weekday_7T == days_chinese_english$chinese[i] & !is.na(demographics$weekday_7T), 'weekday_7T'] <- days_chinese_english$english[i]
}

# Convert scan times that include the date to time of day
scan_times_3T <- gsub(".* (.+)", "\\1", as.character(spreadsheet$`3T MRI Scanning Time`))
scan_times_7T <- gsub(".* (.+)", "\\1", as.character(spreadsheet$`7T MRI Scanning Time`))

# Add to demographics
demographics$scan_times_3T <- scan_times_3T
demographics$scan_times_7T <- scan_times_7T

# Add boolean for whether scanning took place
demographics$did_3T_scan <- !is.na(demographics$weekday_3T)
demographics$did_7T_scan <- !is.na(demographics$weekday_7T)

# /* 
# ----------------------------- Questionnaire data ---------------------------
# */
# Read in the data from the subjects that had complete 2 forgotten questions on paper
# This is strange all of the sudden the code below doesn't work any more
forgottenQuestion_data <- read_excel(paste0(path2data, "Questionnaire.xlsx"))

# Current task folder
taskFolder <- paste0(path2data, "QP_preScanner/")

# Get the input .csv which is used to display the questions so it can be analysed easier. 
inputData_questions  <- read.table(path2input_questions, sep = "\t", header = TRUE, quote = "", encoding = "UTF-8", stringsAsFactors = FALSE)
# The question below was at first forgotten and was later added
forgottenQuestion_ids       <- c("turnbased2_last_year", "turnbased2_before_last_year")
inputData_questions_reduced <- inputData_questions[!(inputData_questions$questionID %in% forgottenQuestion_ids), ]

# Load the subject IDs for this analysis
subjIDs    <- list.files(taskFolder)
allFiles   <- paste0(taskFolder, subjIDs, '/S003/trial_results.csv')
question_data <- do.call(rbind, lapply(allFiles, read.csv, quote = "", encoding = "UTF-8", stringsAsFactors = FALSE))

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$HCnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  question_data$ppid[question_data$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$HCnum == subjIDs[i]]
}

# Relabel to subject & remove ppid & controller_...
question_data$subject <- question_data$ppid
question_data$ppid    <- NULL 
question_data$controller_mouse_screen_location_0 <- NULL
question_data$controller_mouse_screen_location_1 <- NULL

### Add the information form the input file to the question_data
# For this looping through subject variable and based on whether they completed 
# the full questionnaire (53 items) or the one where we forgot two questions (51 items)
question_data$source     <- NA
question_data$questionID <- NA
question_data$english    <- NA

# Subjects to loop over
subjects <- unique(question_data$subject)

# Loop
for(i in 1:length(subjects)){
  if(sum(question_data$subject == subjects[i]) == 51){
    # Case where two questions for not implemented yet
    question_data[question_data$subject == subjects[i], 'source']     <- inputData_questions_reduced$source
    question_data[question_data$subject == subjects[i], 'questionID'] <- inputData_questions_reduced$questionID
    question_data[question_data$subject == subjects[i], 'english']    <- inputData_questions_reduced$english
    
  } else if(sum(question_data$subject == subjects[i]) == 53){
    # Full questionnaire
    question_data[question_data$subject == subjects[i], 'source']     <- inputData_questions$source
    question_data[question_data$subject == subjects[i], 'questionID'] <- inputData_questions$questionID
    question_data[question_data$subject == subjects[i], 'english']    <- inputData_questions$english
  } else
    stop(paste0("Something wrong happened as there are ",  sum(question_data$subject == subjects[i]), " items for subject ", subjects[i]))
}

# Get the question info
inputData_questions    <- inputData_questions[!is.na(inputData_questions$questionID ), ] 
question1 <- inputData_questions[inputData_questions$questionID == forgottenQuestion_ids[1], 'question']
question2 <- inputData_questions[inputData_questions$questionID == forgottenQuestion_ids[2], 'question']
english1  <- inputData_questions[inputData_questions$questionID == forgottenQuestion_ids[1], 'english']
english2  <- inputData_questions[inputData_questions$questionID == forgottenQuestion_ids[2], 'english']
sNA <- c(NA, NA)

# Loop to add the extra questions
for(i in 1:nrow(forgottenQuestion_data)){
  # Get some of the values that are too long to fit below
  vals <- c(forgottenQuestion_data$`Turns-strategy_during_the_past_year`[i], forgottenQuestion_data$`Turns-strategy_before_the_past_year`[i])
  subj <- lookupTable$anonKey[lookupTable$Rnum == forgottenQuestion_data$R_number[i]]
  
  # Add to tempDF
  tempDF <- data.frame(experiment = c("QP_preScanner", "QP_preScanner"),
                       session_num = c(3, 3),
                       trial_num = sNA,
                       block_num = sNA,
                       trial_num_in_block = sNA,
                       start_time = sNA,
                       end_time = sNA,
                       trialType = sNA,
                       question = c(question1, question2),
                       options = sNA,
                       minimumDuration = sNA,
                       value = vals,
                       subject = c(subj, subj),
                       source = c("Video Game Questionnaire", "Video Game Questionnaire"), 
                       questionID = c(forgottenQuestion_ids[1], forgottenQuestion_ids[2]),
                       english = c(english1, english2))
  
  # Add tempDF to question_data
  question_data <- rbind(question_data, tempDF)
}

# Order by subject to make it more sensible
question_data <- question_data[order(question_data$subject), ]

# /* 
# ----------------------------- Diamond World ---------------------------
# */
# Current task folder
taskFolder <- paste0(path2data, "DiamondWorld_standard/")

### DW trial_results
# Load the subject IDs for this analysis
subjIDs    <- list.files(taskFolder)
allFiles   <- paste0(taskFolder, subjIDs, '/S003/trial_results.csv')
DW_trial_results <- do.call(rbind, lapply(allFiles, read.csv, quote = ""))

# Remove columns
DW_trial_results$first_person_controller_movement_location_0 <- NULL
DW_trial_results$first_person_controller_movement_location_1 <- NULL
DW_trial_results$raytracker_ObjectsOnScreenTracker_location_0 <- NULL
DW_trial_results$raytracker_ObjectsOnScreenTracker_location_1 <- NULL

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$HCnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  DW_trial_results$ppid[DW_trial_results$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$HCnum == subjIDs[i]]
}

# Relabel to subject & remove ppid
DW_trial_results$subject <- DW_trial_results$ppid
DW_trial_results$ppid    <- NULL 

### DW log entries
# Create empty list
tempList <- list()

# Loop through all subject folders to load each log file
for(subj in 1:length(subjIDs)){
  # Not this code will only work correctly if the folder only contains position tracker
  # Get all files
  logFile  <- paste0(taskFolder, subjIDs[subj], '/S003/session_info/log.csv') 
  
  # Get list of DFs
  tempLog <- lapply(logFile, read.csv)
  
  # Add index as subject and then bind together
  tempList[[subj]] <- as.data.frame(tempLog)
  
}

# Add index as subject and then bind together into one data frame
DW_logEntries <- custom_rbinder_addIndexColumn(tempList, "subject")

# Add ppid to the data frame
## Get subject rle
subject_rle <- rle(DW_logEntries$subject)

# Add ppid to data frame
DW_logEntries$ppid <- rep(subjIDs, times = subject_rle$lengths) # Repeat each R number in subjIDs according to the run length of subject
# This works because the subjects are loaded in the order of subject_rle

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  DW_logEntries$ppid[DW_logEntries$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$HCnum == subjIDs[i]]
}

########## DW position data
# Load the subject IDs for this analysis
subjIDs    <- list.files(taskFolder)

# Loop through all subjects
# Create empty list
tempList <- list()

# Loop through all subject folders
for(i in 1:length(subjIDs)){
  # Not this code will only work correctly if the folder only contains position tracker
  # Get all files
  path2tracker <- paste0(taskFolder, subjIDs[i], '/S003/trackers/') 
  allTrackers  <- paste0(path2tracker, 'first_person_controller_movement_T001.csv')
  
  # Get list of DFs
  tempTracker <- lapply(allTrackers, read.csv)
  
  # Add index as subject and then bind together
  tempList[[i]] <- custom_rbinder_addIndexColumn(tempTracker, "trial")
}

# Bind this list to 1 DF with subject index
DW_position_data <- custom_rbinder_addIndexColumn(tempList, "subject")

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$HCnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all subject and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  DW_position_data$subject[DW_position_data$subject == which(subjIDs == subjIDs[i])] <- lookupTable$anonKey[lookupTable$HCnum == subjIDs[i]]
}

# Down sample so we can actually analyse the data
DW_position_data <- DW_position_data[seq(from = 1, to = nrow(DW_position_data), by = 10), ]

# /* 
# ----------------------------- OLM: 3T session (practice) ---------------------------
# */
### OLM: 3T session (practice) trial_results
# Current task folder
taskFolder <- paste0(path2data, "OLM_practice/")

# Load the subject IDs for this analysis
subjIDs    <- list.files(taskFolder)
allFiles   <- paste0(taskFolder, subjIDs, '/S003/trial_results.csv')
OLM_3T_practice_trial_results <- do.call(rbind, lapply(allFiles, read.csv, quote = ""))

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$HCnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  OLM_3T_practice_trial_results$ppid[OLM_3T_practice_trial_results$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$HCnum == subjIDs[i]]
}

# Relabel to subject & remove ppid & player_movement_location_0....
OLM_3T_practice_trial_results$subject <- OLM_3T_practice_trial_results$ppid
OLM_3T_practice_trial_results$ppid    <- NULL 
OLM_3T_practice_trial_results$player_movement_location_0 <- NULL
OLM_3T_practice_trial_results$player_movement_location_1 <- NULL

# /* 
# ----------------------------- OLM: 3T session ---------------------------
# */
### OLM: 3T session trial_results
# Current task folder
taskFolder <- paste0(path2data, "OLM_grassy/")

# Load the subject IDs for this analysis
subjIDs    <- list.files(taskFolder)
subjIDs    <- subjIDs[subjIDs %in% lookupTable$Rnum[!(lookupTable$anonKey %in% c("00MCAJ", "U2TQ1G"))]] # This is necessary because two participants were excluded for 3T

#_____________________________
# Needed to be repaired. This is only a temporary solution
subjIDs    <- subjIDs[subjIDs %in% lookupTable$Rnum[!(lookupTable$anonKey %in% c("Z64A18"))]] 
#_____________________________

allFiles   <- paste0(taskFolder, subjIDs, '/S003/trial_results.csv')
OLM_3T_trial_results <- do.call(rbind, lapply(allFiles, read.csv, quote = ""))

# Remove rows with NA values
OLM_3T_trial_results <- na.omit(OLM_3T_trial_results)

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$Rnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  OLM_3T_trial_results$ppid[OLM_3T_trial_results$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$Rnum == subjIDs[i]]
}

# Relabel to subject & remove ppid & player_movement_location_0....
OLM_3T_trial_results$subject <- OLM_3T_trial_results$ppid
OLM_3T_trial_results$ppid    <- NULL
OLM_3T_trial_results$player_movement_location_0 <- NULL
OLM_3T_trial_results$player_movement_location_1 <- NULL

### OLM: 3T session log entries
### DW log entries
# Create empty list
tempList <- list()

# Loop through all subject folders to load each log file
for(subj in 1:length(subjIDs)){
  # Not this code will only work correctly if the folder only contains position tracker
  # Get all files
  logFile  <- paste0(taskFolder, subjIDs[subj], '/S003/session_info/log.csv') 
  
  # Get list of DFs
  tempLog <- lapply(logFile, read.csv)
  
  # Add index as subject and then bind together
  tempList[[subj]] <- as.data.frame(tempLog)
}

# Add index as subject and then bind together into one data frame
OLM_3T_logEntries <- custom_rbinder_addIndexColumn(tempList, "subject")

# Add ppid to the data frame
## Get subject rle
subject_rle <- rle(OLM_3T_logEntries$subject)

# Add ppid to data frame
OLM_3T_logEntries$ppid <- rep(subjIDs, times = subject_rle$lengths) # Repeat each R number in subjIDs according to the run length of subject
# This works because the subjects are loaded in the order of subject_rle

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  OLM_3T_logEntries$ppid[OLM_3T_logEntries$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$Rnum == subjIDs[i]]
}

### OLM: 3T session position data
# Current task folder
taskFolder <- paste0(path2data, "OLM_grassy/")

# Loop through all subjects
# Create empty list
tempList <- list()

# Loop through all subject folders
for(i in 1:length(subjIDs)){
  # Not this code will only work correctly if the folder only contains position tracker
  # Get all files
  path2tracker <- paste0(taskFolder, subjIDs[i], '/S003/trackers/') 
  allTrackers  <- paste0(path2tracker, list.files(path2tracker))
  
  # Get list of DFs
  tempTracker <- lapply(allTrackers, read.csv)
  
  # Add index as subject and then bind together
  tempList[[i]] <- custom_rbinder_addIndexColumn(tempTracker, "trial")
}

# Bind this list to 1 DF with subject index
OLM_3T_position_data <- custom_rbinder_addIndexColumn(tempList, "subject")

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$Rnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all subject and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  OLM_3T_position_data$subject[OLM_3T_position_data$subject == which(subjIDs == subjIDs[i])] <- lookupTable$anonKey[lookupTable$Rnum == subjIDs[i]]
}

# Subset data to moving sections only
OLM_3T_position_data$moving <- ifelse(OLM_3T_position_data$moving == "True", TRUE, FALSE)

# /* 
# ----------------------------- OLM: 7T session (practice) ---------------------------
# */
# ### OLM: 3T session (practice) trial_results
# # Current task folder
# taskFolder <- paste0(path2data, "OLM_practice/")
# 
# # Load the subject IDs for this analysis
# subjIDs    <- list.files(taskFolder)
# 
# # Now check which of these already completed the 7T scan
# completed7T <- c()
# for(i in 1:length(subjIDs)){
#   completed7T[i] <- any(list.files(paste0(taskFolder, subjIDs[i])) == "S007")
# }
# 
# # Subset to those who completed
# subjIDs <- subjIDs[completed7T]
# 
# # Load the data
# allFiles   <- paste0(taskFolder, subjIDs, '/S007/trial_results.csv')
# OLM_7T_practice_trial_results <- do.call(rbind, lapply(allFiles, read.csv, quote = ""))
# 
# # Give error if a subjID cannot be found in the lookUpTable
# if(mean(subjIDs %in% lookupTable$HCnum) != 1){
#   stop("1 & more subjects not found in the look up table.")
# }
# 
# # Loop through all ppid and replace the subject number with the anonymised version
# for(i in 1:length(subjIDs)){
#   OLM_7T_practice_trial_results$ppid[OLM_7T_practice_trial_results$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$HCnum == subjIDs[i]]
# }
# 
# # Relabel to subject & remove ppid & player_movement_location_0....
# OLM_7T_practice_trial_results$subject <- OLM_7T_practice_trial_results$ppid
# OLM_7T_practice_trial_results$ppid    <- NULL 
# OLM_7T_practice_trial_results$player_movement_location_0 <- NULL
# OLM_7T_practice_trial_results$player_movement_location_1 <- NULL

# /* 
# ----------------------------- OLM: 7T session ---------------------------
# */
### OLM: 7T session trial_results
# Current task folder
taskFolder <- paste0(path2data, "OLM_grassy/")

# Load the subject IDs for this analysis
subjIDs    <- list.files(taskFolder)

# Now check which of these already completed the 7T scan
completed7T <- c()
for(i in 1:length(subjIDs)){
  completed7T[i] <- any(list.files(paste0(taskFolder, subjIDs[i])) == "S007")
}

# Subset to those who completed
subjIDs <- subjIDs[completed7T]

# Load the results
allFiles   <- paste0(taskFolder, subjIDs, '/S007/trial_results.csv')
OLM_7T_trial_results <- do.call(rbind, lapply(allFiles, read.csv, quote = ""))

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$Rnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  OLM_7T_trial_results$ppid[OLM_7T_trial_results$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$Rnum == subjIDs[i]]
}

# Relabel to subject & remove ppid & player_movement_location_0....
OLM_7T_trial_results$subject <- OLM_7T_trial_results$ppid
OLM_7T_trial_results$ppid    <- NULL
OLM_7T_trial_results$player_movement_location_0 <- NULL
OLM_7T_trial_results$player_movement_location_1 <- NULL

# Exclude trial that is not completed (e.g. has NA for end_x)
OLM_7T_trial_results <- OLM_7T_trial_results[!is.na(OLM_7T_trial_results$end_x), ]

### OLM: 7T session log entries
### DW log entries
# Create empty list
tempList <- list()

# Loop through all subject folders to load each log file
for(subj in 1:length(subjIDs)){
  # Not this code will only work correctly if the folder only contains position tracker
  # Get all files
  logFile  <- paste0(taskFolder, subjIDs[subj], '/S007/session_info/log.csv') 
  
  # Get list of DFs
  tempLog <- lapply(logFile, read.csv)
  
  # Add index as subject and then bind together
  tempList[[subj]] <- as.data.frame(tempLog)
  
}

# Add index as subject and then bind together into one data frame
OLM_7T_logEntries <- custom_rbinder_addIndexColumn(tempList, "subject")

# Add ppid to the data frame
## Get subject rle
subject_rle <- rle(OLM_7T_logEntries$subject)

# Add ppid to data frame
OLM_7T_logEntries$ppid <- rep(subjIDs, times = subject_rle$lengths) # Repeat each R number in subjIDs according to the run length of subject
# This works because the subjects are loaded in the order of subject_rle

# Loop through all ppid and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  OLM_7T_logEntries$ppid[OLM_7T_logEntries$ppid == subjIDs[i]] <- lookupTable$anonKey[lookupTable$Rnum == subjIDs[i]]
}

### OLM: 7T session position data
# Loop through all subjects
# Create empty list
tempList <- list()

# Loop through all subject folders
for(i in 1:length(subjIDs)){
  # Not this code will only work correctly if the folder only contains position tracker
  # Get all files
  path2tracker <- paste0(taskFolder, subjIDs[i], '/S007/trackers/') 
  allTrackers  <- paste0(path2tracker, list.files(path2tracker))
  
  # Get list of DFs
  tempTracker <- lapply(allTrackers, read.csv)
  
  # Add index as subject and then bind together
  tempList[[i]] <- custom_rbinder_addIndexColumn(tempTracker, "trial")
}

# Bind this list to 1 DF with subject index
OLM_7T_position_data <- custom_rbinder_addIndexColumn(tempList, "subject")

# Give error if a subjID cannot be found in the lookUpTable
if(mean(subjIDs %in% lookupTable$Rnum) != 1){
  stop("1 & more subjects not found in the look up table.")
}

# Loop through all subject and replace the subject number with the anonymised version
for(i in 1:length(subjIDs)){
  OLM_7T_position_data$subject[OLM_7T_position_data$subject == which(subjIDs == subjIDs[i])] <- lookupTable$anonKey[lookupTable$Rnum == subjIDs[i]]
}

# Add trialType to position data
# Get how often the variables that will be added to the movement DF has to be repeated
trial_rle <- rle(OLM_7T_position_data$trial) # Get run length encoding so we know how often we need to repeat the trial type and other variables

# Add the trial type to the position data
OLM_7T_position_data$trialType <- rep(OLM_7T_trial_results$trialType, times = trial_rle$lengths) # Repeat trial type and add to position data

# Subset data to moving sections only
OLM_7T_position_data$moving <- ifelse(OLM_7T_position_data$moving == "True", TRUE, FALSE)

# /* 
# ----------------------------- Overview table ---------------------------
# */
# Create over view table
overview_table <- data.frame(subject = lookupTable$anonKey, Rnum = lookupTable$Rnum)

# Check if data is present & how many trials
# Completed
overview_table$demographics_presented     <- NA
overview_table$DW_completed               <- NA
overview_table$questions_completed        <- NA
overview_table$questions_trials           <- NA
overview_table$OLM_3T_practice_completed  <- NA
overview_table$OLM_3T_practice_trials     <- NA
overview_table$OLM_3T_completed           <- NA
overview_table$OLM_3T_trials              <- NA
overview_table$OLM_7T_completed           <- NA
overview_table$OLM_3T_trials              <- NA
overview_table$OLM_7T_trials              <- NA

# Loop through all subjects
for(i in 1:nrow(overview_table)){
  overview_table[i, 'demographics_presented'] <- any(demographics$subject == overview_table[i, 'subject'])
  
  overview_table[i, 'DW_completed'] <- any(DW_trial_results$subject == overview_table[i, 'subject'])
  
  overview_table[i, 'questions_completed'] <- any(question_data$subject == overview_table[i, 'subject'])
  overview_table[i, 'questions_trials'] <- sum(question_data$subject == overview_table[i, 'subject'])
  
  overview_table[i, 'OLM_3T_practice_completed'] <- any(OLM_3T_practice_trial_results$subject == overview_table[i, 'subject'])
  overview_table[i, 'OLM_3T_practice_trials'] <- sum(OLM_3T_practice_trial_results$subject == overview_table[i, 'subject'])
  
  overview_table[i, 'OLM_3T_completed'] <- any(OLM_3T_trial_results$subject == overview_table[i, 'subject'])
  overview_table[i, 'OLM_3T_trials'] <- sum(OLM_3T_trial_results$subject == overview_table[i, 'subject'])
  
  overview_table[i, 'OLM_7T_completed'] <- any(OLM_7T_trial_results$subject == overview_table[i, 'subject'])
  overview_table[i, 'OLM_7T_trials'] <- sum(OLM_7T_trial_results$subject == overview_table[i, 'subject'])
}

# Sort by Rnum & remove row names
overview_table <- overview_table[order(overview_table$Rnum), ]
row.names(overview_table) <- NULL

# /* 
# ----------------------------- Remove any incomplete data sets ---------------------------
# */
# 3T
exclude <- overview_table$subject[overview_table$OLM_3T_trials != 36]
#OLM_3T_trial_results <- OLM_3T_trial_results[!(OLM_3T_trial_results$subject %in% exclude), ]
#OLM_3T_logEntries    <- OLM_3T_logEntries[!(OLM_3T_logEntries$subject %in% exclude), ]
#OLM_3T_position_data <- OLM_3T_position_data[!(OLM_3T_position_data$subject %in% exclude), ]

# 7T Also include the participant that has 102 trials
exclude <- overview_table$subject[overview_table$OLM_7T_trials != 108 & overview_table$OLM_7T_trials != 102]
#OLM_7T_trial_results <- OLM_7T_trial_results[!(OLM_7T_trial_results$subject %in% exclude), ]
#OLM_7T_logEntries    <- OLM_7T_logEntries[!(OLM_7T_logEntries$subject %in% exclude), ]
#OLM_7T_position_data <- OLM_7T_position_data[!(OLM_7T_position_data$subject %in% exclude), ]

# /* 
# ----------------------------- Save data ---------------------------
# */
# Total exclusion:
## overview table
# Save as .RData
save(list = c("overview_table"), file = paste0(path2output, "overview_table.RData"))

# .txt files
write.table(overview_table, file = paste0(path2output, "overview_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

## Demographics
# Save as .RData
save(list = c("demographics"), file = paste0(path2output, "demographics.RData"))

# .txt files
write.table(demographics, file = paste0(path2output, "demographics.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

## DW
# Save as .RData
save(list = c("DW_trial_results", "DW_logEntries", "DW_position_data"), file = paste0(path2output, "DW_all_data.RData"))

# .txt files
write.table(DW_trial_results, file = paste0(path2output, "DW_trial_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(DW_logEntries,    file = paste0(path2output, "DW_logEntries.txt"),    sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(DW_position_data, file = paste0(path2output, "DW_position_data.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

## Questionnaires
# Save as .RData
save(list = c("question_data"), file = paste0(path2output, "question_data.RData"))

# .txt files
write.table(question_data, file = paste0(path2output, "question_data.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

## OLM 3T practice
# Save as .RData
save(list = c("OLM_3T_practice_trial_results"), file = paste0(path2output, "OLM_3T_practice_trial_results.RData"))

# .txt files
write.table(OLM_3T_practice_trial_results, file = paste0(path2output, "OLM_3T_practice_trial_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

## OLM 3T
# Save as .RData
save(list = c("OLM_3T_trial_results", "OLM_3T_logEntries", "OLM_3T_position_data"), file = paste0(path2output, "OLM_3T_all_data.RData"))

# .txt files
write.table(OLM_3T_trial_results, file = paste0(path2output, "OLM_3T_trial_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(OLM_3T_logEntries,    file = paste0(path2output, "OLM_3T_logEntries.txt"),    sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(OLM_3T_position_data, file = paste0(path2output, "OLM_3T_position_data.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

## OLM 7T practice
###################### ADDD at some point

## OLM 7T
# Save as .RData
save(list = c("OLM_7T_trial_results", "OLM_7T_logEntries", "OLM_7T_position_data"), file = paste0(path2output, "OLM_7T_all_data.RData"))

# .txt files
write.table(OLM_7T_trial_results, file = paste0(path2output, "OLM_7T_trial_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(OLM_7T_logEntries,    file = paste0(path2output, "OLM_7T_logEntries.txt"),    sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(OLM_7T_position_data, file = paste0(path2output, "OLM_7T_position_data.txt"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")