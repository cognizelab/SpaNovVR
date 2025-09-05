#!/bin/bash

# Define where the study and task folder is
export DATA=/home1/Jaquent/Datasets/VP00157_HCP
export TASK=tfMRI_OLMe
export PARAM=hp200_s2
export subject_suffix="_7TSESS1"
export toolboxes_path=/home1/Jaquent/Toolboxes
export analysisName="SpaNov_gradient_6lvl_cue-delay_smo2_MSMAll"
export full_analysis_name="OLMe_7T_${analysisName}"

# Number of contrasts
lvl2_COPES=(1)
num_lvl1_cope=14

# If it does not already exist, create log folder because slurm needs it
mkdir -p log

# Submit each job to SLURM
for lvl1_COPE in `seq 1 ${num_lvl1_cope}`; do
	for lvl2_COPE in ${lvl2_COPES[*]}; do
		export lvl2_COPE
		export lvl1_COPE
		sbatch /home1/Jaquent/scripts/analysis_utility/PALM_slurm/PALM_jobs_multi_2ndLevel_cope_MSMAll.slurm
	done
done

# Check if the jobs are submted
sleep 4
squeue
