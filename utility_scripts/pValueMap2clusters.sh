#!/bin/bash
# This scripts takes a p-Value map and finds the clusters in them by creating a binary mask first.
# Code adapted from Deniz Vatansever

# New line
echo ""
echo "Finding clusters"

# Arguments passed on
pValueMap=$1
echo "P-value for which to find clusters: ${pValueMap}"

# Check if file is supposed to be saved in a different place
if [[ ! -z $2 ]];then
    pThreshold=$2
    echo "Use custom significance cut-off of ${pThreshold}"
else
    pThreshold=1.301
    echo "Use default significance cut-off of ${pThreshold}"
fi

# Create names for binary mask and cluster map
binMask=${pValueMap/.dscalar.nii/_bin.dscalar.nii}
clusterMap=${pValueMap/.dscalar.nii/_clusters.dscalar.nii}
echo "Binary mask will be called: ${binMask}"

# Fixed parameters
path2atlas="/media/alex/work/Seafile/imaging_results/HCP_S1200_Atlas/HCP_S1200_Atlas"

### To create a binarised mask
wb_command -cifti-math "mask > ${pThreshold}" -var mask ${pValueMap} ${binMask}

### To find clusters within the binarised mask
wb_command -cifti-find-clusters ${binMask} 0 0 0 0 COLUMN ${clusterMap} -left-surface "${path2atlas}/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii" -right-surface "${path2atlas}/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii" -merged-volume
