#!/bin/bash
################################################################ 
# This script will run a HCP-style fMRI task analysis on the data. 
# Pre-requisites for running this script: 
# 1) the data & event files were copied into the correct folders.
# 2) The .fsf files are ready. 
# 3) sessions_hcp_batch.txt is edited with correct paths.

# Analysis-specific things that need to be specified
export subject2analyse_file="subjects2analyse.txt"
export data_folder="/home1/Jaquent/Datasets/VP00157_HCP"
export jobname="SpaNovGradient"
export subject_suffix="_7TSESS1"
export task_name="OLMe"
export scanner="7T"
export fullName="SpaNov_gradient_6lvl_smo4_MSMAll"
export fsf_files_folder="/home1/Jaquent/research_projects/OLM_project/analsyis/imaging/fsl_templates/${task_name}_${scanner}_${fullName}"
export event_files="${task_name}_${scanner}_${fullName}"

# Analysis parameter
export smoothing=4
export filtering=200

# SLURM paramteres
export memory="5GB"
export time="2-00:00:00"
export email="jaquent@fudan.edu.cn"

# General things that need to be specified
export path_2_qunex_container_script="/home1/Jaquent/Toolboxes/qunex"
export qunex_container_file="/home1/Jaquent/Toolboxes/my_images/qunex_suite-0.95.2.sif"

# Fixed parameters that should need no change
# Create folder names
# E.g. tfMRI_OLMr_RUN1_AP
export level1_run1_folder="tfMRI_${task_name}_RUN1_AP"
export level1_run2_folder="tfMRI_${task_name}_RUN2_AP"
# E.g. tfMRI_OLMr
export level2_folder="tfMRI_${task_name}"

# Create .fsf file names
# E.g. tfMRI_OLMr_RUN1_AP_Accuracy
export level1_run1_fsf="tfMRI_${task_name}_RUN1_AP_${fullName}"
export level1_run2_fsf="tfMRI_${task_name}_RUN2_AP_${fullName}"
# E.g. tfMRI_OLMr_Accuracy
export level2_fsf="tfMRI_${task_name}_${fullName}"

# Print out the parameter used so mistakes can be spotted
echo "Parameter that have been used. Check carefully..."
echo "fsf_files_folder = ${fsf_files_folder}"
echo "level1_run1_folder = ${level1_run1_folder}"
echo "level2_folder = ${level2_folder}"
echo "level1_run1_fsf = ${level1_run1_fsf}"
echo "level2_fsf = ${level2_fsf}"
echo "event_files = ${event_files}"

################################################################ 
# This section parses and prepare the information specified above
# Read the .txt that specifies the subject that needs to be analysed
subject_string=""
while read -r line; do
  subject_string+="$line\n" # Add each line to the end of the string
done < "$subject2analyse_file"

# Remove the \n from the subject string
subject_string=${subject_string//[\\n]/}

# Split string to array
IFS=',' read -r -a subjects <<< "$subject_string"

# Number of subjects
subject_num=${#subjects[@]}

# Print to console
echo "##################################################"
echo "${subject_num} subjects were specified in the subject file." 

################################################################ 
# This section copies the .fsf & EVs files into correct locations
echo "##################################################"
echo "Copying the .fsf & EV files into the correct locations." 

# Delete all .fsf files
# Source https://www.cyberciti.biz/faq/linux-unix-how-to-find-and-remove-files/
#find ${data_folder} -type f -name "*.fsf" -exec rm -vf {} \;

# load fsl & miniconda3
module load fsl
module load miniconda3

# index
i=1

# Loop over subjects
for SUBJ in "${subjects[@]}"; do
	# Loop over runs
	for BOLD in ${level1_run1_folder} ${level1_run2_folder}; do
		# Added because I get cp: cannot create regular file message while running it the first time.
		mkdir -p ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${BOLD}

		# Removing old files to clean up
		rm -rf ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${BOLD}/EVs

		# Copying the level 1 .fsf file of this run
		cp ${fsf_files_folder}/${BOLD}_${fullName}_hp200_s4_level1.fsf ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${BOLD}/

		# Replacing the placeholder volumes (123) with the correct number of volumes
		NVOLS=$(fslnvols ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${BOLD}/${BOLD}.nii.gz)
		sed -z -i "s|123|${NVOLS}|1" ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${BOLD}/${BOLD}_${fullName}_hp200_s4_level1.fsf
	done

	# Copying the EVs into the correct folder
	cp -R ${data_folder}/sessions/inbox/events/${event_files}/sub-${SUBJ}/ses-01/run-01/EVs ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${level1_run1_folder}/
	cp -R ${data_folder}/sessions/inbox/events/${event_files}/sub-${SUBJ}/ses-01/run-02/EVs ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${level1_run2_folder}/

	# Make dir again
	mkdir -p ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${level2_folder}

	# Copying level 2 .fsf file
	cp -R ${fsf_files_folder}/${level2_folder}/${level2_fsf}_hp200_s4_level2.fsf ${data_folder}/sessions/${SUBJ}${subject_suffix}/hcp/${SUBJ}${subject_suffix}/MNINonLinear/Results/${level2_folder}/${level2_fsf}_hp200_s4_level2.fsf

	# Print status to console
	echo "Processing of ${SUBJ} completed. Progress: ${i}/${subject_num} subjects."

	# Updating index
	((i+=1))

done

echo "Copying completed."

################################################################
# Create the sessionids based on subjects
separator=","
sessionids=""

# Loop through each element in the subjects array and append it to the output string
for (( i=0; i<${#subjects[@]}; i++ )) ; do
    sessionids="${sessionids}${subjects[$i]}${subject_suffix}${separator}"
done

# Remove the last separator character from the output string
# Not sure how this works but it does work.
export sessionids=${sessionids%$separator}

# Print to console
echo "##################################################"
echo "Submitting qunex job." 

# Trying to submit the following particpants:
echo "Submitting: ${sessionids}"

# This section runs the qunex command.
python3 $path_2_qunex_container_script/qunex_container hcp_task_fmri_analysis \
--sessionsfolder="${data_folder}/sessions/" \
--sessions="${data_folder}/processing/sessions_hcp_batch.txt" \
--sessionids="${sessionids}" \
--overwrite="yes" \
--bash_pre="module load singularity" \
--bind="/public/home2/VNLab/" \
--hcp_task_lvl1tasks="${level1_run1_folder}@${level1_run2_folder}" \
--hcp_task_lvl2task="${level2_folder}" \
--hcp_task_lvl1fsfs="${level1_run1_fsf}@${level1_run2_fsf}" \
--hcp_task_lvl2fsf="${level2_fsf}" \
--hcp_task_highpass="${filtering}" \
--hcp_bold_final_smoothFWHM="${smoothing}" \
--hcp_task_confound="Confound_Parameters.txt" \
--hcp_task_procstring="MSMAll_hp0_clean" \
--container=$qunex_container_file \
--scheduler="SLURM,mem-per-cpu=${memory},time=${time}, partition=blade, jobname=${jobname},mail-user=${email}, mail-type=ALL, exclude=bnode05"

# For parcellated
#--hcp_task_parcellation="MMP1.0" \
#--hcp_task_parcellation_file="/public/home2/VNLab/Atlases/MMP1.0_210V_Parcellation/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii" \

# Check if jobs are actually submitted
sleep 1
squeue
