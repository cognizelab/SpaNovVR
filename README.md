# Repository for Graded encoding of spatial novelty scales in the human brain
Authors: Jörn Alexander Quent, Liangyue Song, Xinyu Liang, Yueting Su, Wenwen Yu, He Wang and Deniz Vatansever

To be published in Nature Communications

This repository contains the code for our manuscript investigating gradients of spatial novelty and familiarity in the human hippocampus and beyond.

## Abstract
Successful navigation relies on the ability to process and encode detailed information about our dynamic environments. Beyond familiarity, emerging studies now highlight the crucial role of novelty detection in this process, the precise neural mechanism of which remains poorly understood. Using ultra-high field 7T fMRI, we investigated how the human brain encodes spatial novelty during virtual navigation, with a particular focus on graded representations that follow systematic transitions between novel and familiar spaces. Our results revealed novelty and familiarity specific neural responses within the posterior and anterior poles of the bilateral hippocampus, respectively. On the cortical surface, two separable streams of activity patterns were observed in which regions within the visual and frontoparietal networks showed novelty-specific activity, while somatomotor, ventral attention and default mode regions preferred spatial familiarity. Importantly, we identified a distinct gradient along the hippocampal long axis and demonstrated the extended contribution of the posterior medial cortex to the encoding of spatial novelty scales that were intrinsically coupled with the hippocampal gradient. These findings advance our understanding of how the human brain encodes and processes spatial information, suggesting that graded representations of spatial novelty may serve as a fundamental organizational principle for spatial cognition in the human brain. 

## Keywords
Spatial novelty, hippocampus, gradients, virtual navigation, posterior medial cortex, 7T fMRI

## Additional resources
- [OLM](https://github.com/JAQuent/Object-Location-Memory-Task) the source code & builds to the Object Location Memory task used for this study. Note that we used **Version 2.1.1** of this task for the data collection.
- [Diamond World](https://github.com/JAQuent/DiamondWorld), a foraging task we used to practice movement.
- [Questionnaire Presenter](https://github.com/JAQuent/Questionnaire-Presenter) questionnaire data was collected using this Unity-based program.


## Guide
Here is a detailed guide to the code & files:

- ***3D_hist_blender*** Blender files to create the 3D histograms.
- ***brms_scripts*** folder of all the scripts for the Bayesian hierarchical models. 
	- ***SpaNov_paper_brms_connectivity2.R*** script for running Bayesian hierarchical for the connectivity analysis.
	- ***SpaNov_paper_brms_navTime.R*** script for running Bayesian hierarchical anlaysing navigation times as a function of the number of times an object was presented.
	- ***SpaNov_paper_brms_noveltyScore.R*** script to run the `brms` model of aggregate spatial novelty as a function of time.
	- ***SpaNov_paper_brms_timeSince.R*** script to run the `brms` model of time elapsed between visits to sector as a function of time.
	- ***SpaNov_paper_brms_visits_poisson.R*** script to run the `brms` model of numbe of visits to sector as a function of time modelled with a Poisson distribution.
- ***preparation*** notebook to generate trial sequence.
- ***resting_state_analysis***
	- ***analyse_conn_via_SpaNov_gradient.R*** script to calculate average connectivity between gradients in hippocampus and PMC (+ .slurm script).
	- ***Get_HC_2_cortex_conn_matrices.R*** script to calculate the dense connectivity between each hippocampal voxel and each cortical vertex (+ .slurm script).
- ***tables*** folder that contains the spreadsheets for the supplementary tables. 
- ***tSNR_GS*** scripts to calculate temporal signal-to-noise ratio and global signal for oru data.
- ***utility_scripts***
	- ***anonymiser_script.R*** script to anonymise the participant IDs.
	- ***data_combination_script.R*** script to concatenate and prepare the raw data for analysis.
	- ***hippocampus_normal_point_to_curve_projection.R*** script to calculate how anterior each voxel in the hippocampus is. 
	- ***pValueMap2clusters.sh*** script to find significant clusters (see below for code example).
- ***get_noveltyGradient_slope.R*** script to calculate hippocampal gradient for each participant separately.
- ***OLMe_7T_encoding1_cue-delay_smo4_MSMAll_HCP.sh*** &  ***SOLMe_7T_encoding1_cue-delay_smo4_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which received further smoothing (4 mm), was used to the GLM analyses for encoding success.
- ***SpaNov_cortical_gradient_analysis_min.R*** scripts to run the depth-first algorithm on the cortical gradient maps.
- ***SpaNov_event_file_creation_encoding.Rmd*** (+ .html) the script to create the event files for the encoding success analysis.
- ***SpaNov_event_file_creation_gradient.Rmd*** (+ .html) the script to create the event files for the spatial novelty analysis.
- ***SpaNov_event_file_creation_gradient_dist2centre.Rmd*** (+ .html) the script to create the event files for the spatial novelty analysis, where each events' novelty score is corrected for centrality.
- ***SpaNov_gradient_6lvl_cue-delay_smo2_MSMAll_HCP.sh*** &  ***SpaNov_gradient_6lvl_cue-delay_smo2_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which did not receive further smoothing (2 mm), was used to the gradient analyses for spatial novelty.
- ***SpaNov_gradient_6lvl_cue-delay_smo4_MSMAll_HCP.sh*** &  ***SpaNov_gradient_6lvl_cue-delay_smo4_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which received further smoothing (4 mm), was used to the GLM analysesfor spatial novelty.
- ***SpaNov_gradient_dist2centre-corrected_6lvl_smo4_MSMAll_HCP.sh*** &  ***SpaNov_gradient_dist2centre-corrected_6lvl_smo4_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which received further smoothing (4 mm), was used to the GLM analysesfor spatial novelty, where each events' novelty score is corrected for centrality.
- ***SpaNov_paper_calculate_gradients.R*** script to create the gradient maps.
- ***SpaNov_results.Rmd*** (+ .html) the main script for the analyses.
- ***spin_test_script.py*** & ***spin_test.sh*** scripts to run the spin test to compare the cortical maps.

Example usage of `pValueMap2clusters.sh`:

```
pValueMap="/media/alex/work/Seafile/imaging_results/SpaNov/OLMe_7T_SpaNov_gradient_6lvl_cue-delay_smo4_MSMAll/cope7.feat/stats/vwc/results_lvl2cope1_dat_ztstat_cfdrp_c1.dscalar.nii"
bash /media/alex/work/research_projects/OLM_project/analysis/utility_scripts/pValueMap2clusters.sh ${pValueMap}
```

