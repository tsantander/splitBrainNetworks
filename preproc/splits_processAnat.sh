#!/bin/bash
#-----------------------------------------------------------------------------------------
# Process anatomical data using ANTS:
#	1) Spatially-adaptive denoising
#	2) antsCorticalThickness for high-quality segmentation / brain extraction
#	3) Post-process segmentation posteriors 
#_________________________________________________________________________________________
#
# Part of the preprocessing and analysis code for:
# Santander et al. (2025). Full interhemispheric integration sustained by a fraction of 
# posterior callosal fibers. Proc Natl Acad Sci USA, 122(43) e2520190122. 
# https://doi.org/10.1073/pnas.2520190122
#
# Github repo: https://github.com/tsantander/splitBrainNetworks
#_________________________________________________________________________________________

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

# Preliminary setup.
#-----------------------------------------------------------------------------------------
# Patient ID should be passed as an input argument, so grab that.

patID=$1

# Make sure we're in the right directory.

cd "/zwork/tyler/splitMRI/data/${patID}"

# Denoise T1.
#-----------------------------------------------------------------------------------------
DenoiseImage -d 3 \
	-i "${PWD}/anat/${patID}_T1w.nii.gz" \
	-o "${PWD}/anat/${patID}_T1w_denoised.nii.gz" \
	-v 1
	
# Run cortical thickness pipeline.
#-----------------------------------------------------------------------------------------
# First specify directory for reference template.

templateDir="/zwork/tyler/splitMRI/template/refTemp"

# Technically we really only need to run the first three stages of antsCorticalThickness 
# in order to get what we need here, but might as well execute the whole pipeline so we 
# have these data for later.

antsCorticalThickness.sh -d 3 \
	-a "${PWD}/anat/${patID}_T1w_denoised.nii.gz" \
	-e "${templateDir}/T_template0.nii.gz" \
	-t "${templateDir}/T_template0_BrainCerebellum.nii.gz" \
	-m "${templateDir}/T_template0_BrainCerebellumProbabilityMask.nii.gz" \
	-f "${templateDir}/T_template0_BrainCerebellumExtractionMask.nii.gz" \
	-p "${templateDir}/Priors/priors%d.nii.gz" \
	-n 6 \
	-x 6 \
	-o "${PWD}/anat/${patID}_T1w_"

# Binarize CSF and WM posteriors for nuisance signal estimation later.
#-----------------------------------------------------------------------------------------
# First CSF.

ThresholdImage 3 "${PWD}/anat/${patID}_T1w_BrainSegmentationPosteriors1.nii.gz" "${PWD}/anat/csfMask.nii.gz" 0.99 1 1 0

# Now WM.

ThresholdImage 3 "${PWD}/anat/${patID}_T1w_BrainSegmentationPosteriors3.nii.gz" "${PWD}/anat/wmMask.nii.gz" 0.99 1 1 0

# Smooth posteriors.
#-----------------------------------------------------------------------------------------
for iSeg in {1..6}; do

    SmoothImage 3 "${PWD}/anat/${patID}_T1w_BrainSegmentationPosteriors${iSeg}.nii.gz" 2 "${PWD}/anat/${patID}_T1w_smoothedPosteriors${iSeg}.nii.gz" 1
    
done

exit 0
