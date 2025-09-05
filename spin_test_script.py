# -*- coding: utf-8 -*-
"""
Running spin test for the spatial novelty paper
================================================
"""

###############################################################################
# Path to the CIFTI files
map1_path = '/media/alex/work/Seafile/imaging_results/SpaNov/OLMe_7T_SpaNov_gradient_6lvl_smo4_MSMAll/cope7.feat/stats/vwc/results_lvl2cope1_dat_ztstat_c1.dscalar.nii'
map2_path = '/media/alex/work/Seafile/imaging_results/SpaNov/OLMe_7T_SpaNov_gradient_6lvl_cue-delay_smo4_MSMAll/cope7.feat/stats/vwc/results_lvl2cope1_dat_ztstat_c1.dscalar.nii'
map3_path = '/media/alex/work/Seafile/imaging_results/SpaNov/OLMe_7T_SpaNov_gradient_dist2centre-corrected_6lvl_smo4_MSMAll/cope7.feat/stats/vwc/results_lvl2cope1_dat_ztstat_c1.dscalar.nii'

# Parameteres
n_perm = 10000

# Libraries
from neuromaps import stats
from neuromaps import nulls
import nibabel as nib
import numpy as np

###############################################################################
def load_cifti_32k_hemispheres(cifti_path):
    """
    Load a CIFTI file and return left and right hemisphere data
    with 32,492 vertices each, preserving medial wall values as NaN.
    """
    print("\n===============================")
    print(f"Loading CIFTI file: {cifti_path}")
    
    # Load the CIFTI file
    img = nib.load(cifti_path)
    data = img.get_fdata()
    print(f"Original data shape: {data.shape}")
    
    # Initialize arrays for both hemispheres with NaN values (medial wall)
    left_hemisphere = np.full(32492, np.nan)
    right_hemisphere = np.full(32492, np.nan)
    
    # Get brain model axis to map data to vertices
    brain_models = img.header.get_axis(1)
    
    # Map CIFTI data to appropriate vertices
    for name, indices, brain_model in brain_models.iter_structures():
        if name == 'CIFTI_STRUCTURE_CORTEX_LEFT':
            left_hemisphere[brain_model.vertex] = data[0, indices]
            # Calculate count properly for slice objects
            if isinstance(indices, slice):
                start, stop, step = indices.indices(data.shape[1])
                count = len(range(start, stop, step or 1))
            else:
                count = len(indices)
            print(f"Mapped {count} values to left hemisphere")
        elif name == 'CIFTI_STRUCTURE_CORTEX_RIGHT':
            right_hemisphere[brain_model.vertex] = data[0, indices]
            # Calculate count properly for slice objects
            if isinstance(indices, slice):
                start, stop, step = indices.indices(data.shape[1])
                count = len(range(start, stop, step or 1))
            else:
                count = len(indices)
            print(f"Mapped {count} values to right hemisphere")
    
    # Verify both hemispheres have correct number of vertices
    assert left_hemisphere.size == 32492, "Left hemisphere vertex count mismatch"
    assert right_hemisphere.size == 32492, "Right hemisphere vertex count mismatch"
    
    # Concatenate them
    brain_map = np.concatenate([left_hemisphere, right_hemisphere])

    print("Successfully loaded and processed CIFTI file")
    print(f"Left hemisphere shape: {left_hemisphere.shape}")
    print(f"Right hemisphere shape: {right_hemisphere.shape}")
    print(f"Concatenated brain map shape: {brain_map.shape}")
    
    return brain_map

###############################################################################
# Extract cortical data from all maps
print("Extracting cortical data from CIFTI files...")
map1_data = load_cifti_32k_hemispheres(map1_path)
map2_data = load_cifti_32k_hemispheres(map2_path)
map3_data = load_cifti_32k_hemispheres(map3_path)

################################################################################
print("\n===============================")
print("Rotating null distribution...")
rotated = nulls.alexander_bloch(map2_data, atlas='fsLR', density='32k', 
                                n_perm=n_perm, seed=1225)
print(rotated.shape)

################################################################################
print("\n===============================")
corr1, pval1 = stats.compare_images(map2_data, map1_data, nulls=rotated)
print(f'Correlation between Map 2 and Map 1: r = {corr1:.02f}, p = {pval1:.04f}')
corr2, pval2 = stats.compare_images(map2_data, map3_data, nulls=rotated)
print(f'Correlation between Map 2 and Map 3: r = {corr2:.02f}, p = {pval2:.04f}')
#Correlation between Map 2 and Map 1: r = 0.98, p = 0.0001
#Correlation between Map 2 and Map 3: r = 0.91, p = 0.0001
