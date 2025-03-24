# Repository for SpaNovVR
Collaborators: Liangyue Song, Yueting Su, Wenwen Yu, He Wang, Xinyu Liang, Deniz Vatansever

Status: submitted

This repository contains the code our manuscript investigating gradients of spatial novelty and familiarity in the human hippocampus and beyond.

## Additional resources:
- [OLM](https://github.com/JAQuent/Object-Location-Memory-Task) the source code & builds to the Object Location Memory task used for this study. Note that we used **Version 2.1.1** of this task for the data collection.
- [Diamond World](https://github.com/JAQuent/DiamondWorld), a foraging task we used to practice movement.
- [Questionnaire Presenter](https://github.com/JAQuent/Questionnaire-Presenter) questionnaire data was collected using this Unity-based program.


## Guide
Here is a detailed guide to the code & files:

- ***anonymiser_script.R*** the script to anonumise the participant IDs.
- ***data_combination_script.R*** the script to concatenate and prepare the data for further analysis.
- ***SpaNov_event_file_creation_gradient.html*** & ***SpaNov_event_file_creation_gradient.html*** the script to create the event files for the analysis.
- ***SpaNov_gradient_6lvl_smo2_MSMAll_HCP.sh*** &  ***SpaNov_gradient_6lvl_smo2_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which did not receive further smoothing (2 mm) was used to the gradient analyses.
- ***SpaNov_gradient_6lvl_smo4_MSMAll_HCP.sh*** &  ***SpaNov_gradient_6lvl_smo4_MSMAll_PALM.sh*** the script to run the HCP-style task fMRI analysis using a QuNex container and the script to run the corresponding the PALM analysis. This version, which received further smoothing (4 mm) was used to the GLM analyses.