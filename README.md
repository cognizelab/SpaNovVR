# Repository for Graded encoding of spatial novelty scales in the human brain
Collaborators: Liangyue Song, Yueting Su, Wenwen Yu, He Wang, Xinyu Liang, Deniz Vatansever

Status: submitted

This repository contains the code for our manuscript investigating gradients of spatial novelty and familiarity in the human hippocampus and beyond.

## Additional resources
- [OLM](https://github.com/JAQuent/Object-Location-Memory-Task) the source code & builds to the Object Location Memory task used for this study. Note that we used **Version 2.1.1** of this task for the data collection.
- [Diamond World](https://github.com/JAQuent/DiamondWorld), a foraging task we used to practice movement.
- [Questionnaire Presenter](https://github.com/JAQuent/Questionnaire-Presenter) questionnaire data was collected using this Unity-based program.


## Guide
Here is a detailed guide to the code & files:

- ***3D_hist_blender*** Blender files to create the 3D histograms.
- ***resting_state_analysis***
	- ***analyse_conn_via_SpaNov_gradient.R*** script to calculate average connectivity between gradients in hippocampus and PMC (+ .slurm script).
	- ***Get_HC_2_cortex_conn_matrices.R*** script to calculate the dense connectivity between each hippocampal voxel and each cortical vertex (+ .slurm script).
- ***tables*** folder that contains the spreadsheets for the supplementary tables. 
- ***utility_scripts***
	- ***anonymiser_script.R*** script to anonymise the participant IDs.
	- ***data_combination_script.R*** script to concatenate and prepare the raw data for analysis.
	- ***hippocampus_normal_point_to_curve_projection.R*** script to calculate how anterior each voxel in the hippocampus is. 
	- ***pValueMap2clusters.sh*** script to find significant clusters (see below for code example).
- ***SpaNov_cortical_gradient_analysis_max.R*** & ***SpaNov_cortical_gradient_analysis_min.R*** scripts to run the depth-first algorithm on the cortical gradient maps.
- ***SpaNov_event_file_creation_gradient.html*** & ***SpaNov_event_file_creation_gradient.html*** the script to create the event files.
- ***SpaNov_gradient_6lvl_smo2_MSMAll_HCP.sh*** &  ***SpaNov_gradient_6lvl_smo2_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which did not receive further smoothing (2 mm), was used to the gradient analyses.
- ***SpaNov_gradient_6lvl_smo4_MSMAll_HCP.sh*** &  ***SpaNov_gradient_6lvl_smo4_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which received further smoothing (4 mm), was used to the GLM analyses.
- ***SpaNov_paper_brms_noveltyScore.R*** script to run the `brms` model of spatial novelty as a function of time.
- ***SpaNov_paper_calculate_gradients.R*** script to create the gradient maps.


Example usage of `pValueMap2clusters.sh`:

```
pValueMap="/media/alex/work/Seafile/imaging_results/SpaNov/OLMe_7T_SpaNov_gradient_6lvl_smo4_MSMAll/cope7.feat/stats/vwc/results_lvl2cope1_dat_ztstat_cfdrp_c1.dscalar.nii"
bash /media/alex/work/research_projects/OLM_project/analysis/utility_scripts/pValueMap2clusters.sh ${pValueMap}
```

