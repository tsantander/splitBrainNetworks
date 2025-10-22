#!/bin/bash
#-----------------------------------------------------------------------------------------
# Spatial registration:
# 	1) EPI --> T1
#	2) T1 --> MNI (using predefined cost-function mask in target space)
#	3) Concatenate EPI --> T1 --> MNI
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

# Patient ID should be passed as an input argument, along with the number of resting-state
# runs and the corresponding CFM, so grab all that.

patID=$1
costFunMask=$2
nRest=$3

# Make sure we're in the right directory.

cd "/zwork/tyler/splitMRI/data/${subID}"

# Register mean EPIs to skull-stripped T1.
#-----------------------------------------------------------------------------------------
pathToT1="/zwork/tyler/splitMRI/data/${patID}/anat"
pathToEPI="/zwork/tyler/splitMRI/data/${patID}/func"

# Grab the skull-stripped T1 and the brain mask.

ssT1="${pathToT1}/${patID}_T1w_ExtractedBrain0N4.nii.gz"

brainMask="${pathToT1}/${patID}_T1w_BrainExtractionMask.nii.gz"

# Dilate extraction mask for coregistration.

ImageMath 3 "${pathToT1}/${patID}_T1w_BrainExtractionMask_dil.nii.gz" MD ${brainMask} 10

coregMaskT1="${pathToT1}/${patID}_T1w_BrainExtractionMask_dil.nii.gz"

if [ $nRest -eq 1 ]; then

	# Bias correct mean EPI (less-aggressive N3).
		
	N3BiasFieldCorrection 3 \
		"${PWD}/meanEPI/meanRestEPI.nii" \
		"${PWD}/meanEPI/meanRestN3.nii.gz"

	meanRestEPI="${PWD}/meanEPI/meanRestN3.nii.gz"
	
	# Quick brain mask estimation.
	
	mri_synthstrip -i ${meanRestEPI} -m "${PWD}/meanEPI/meanRestN3_mask.nii.gz"
	ImageMath 3 "${PWD}/meanEPI/meanRestN3_mask_dil.nii.gz" MD "${PWD}/meanEPI/meanRestN3_mask.nii.gz" 10
	
	coregMaskEPI="${PWD}/meanEPI/meanRestN3_mask_dil.nii.gz"
	
	# Estimate coregistration.
	
	antsRegistration -d 3 \
		-o "${PWD}/meanEPI/mean2anat_rest_" \
		-r [${ssT1},${meanRestEPI},1] \
		-t Translation[0.10] \
		-m Mattes[${ssT1},${meanRestEPI},1,32,Regular,0.30] \
		-c [12000x12000x11110,1e-8,20] \
		-s 4x2x1vox \
		-f 6x4x2 \
		-u 0 \
		-t Rigid[0.10] \
		-m Mattes[${ssT1},${meanRestEPI},1,32,Regular,0.30] \
		-c [12000x12000x11110,1e-8,20] \
		-s 4x2x1vox \
		-f 3x2x1 \
		-u 0 \
		-t Affine[0.10] \
		-m Mattes[${ssT1},${meanRestEPI},1,32,Regular,0.30] \
		-c [12000x12000x11110,1e-8,20] \
		-s 4x2x1vox \
		-f 3x2x1 \
		-u 0 \
		-t SyN[0.20,3,0] \
		-m Mattes[${ssT1},${meanRestEPI},0.5,32] \
		-m CC[${ssT1},${meanRestEPI},0.5,4] \
		-c [100x100x50,-0.01,5] \
		-s 1x0.5x0vox \
		-f 4x2x1 \
		-u 0 \
		-w [0.025,0.975] \
		-z 1 \
		-v 1

				
else
	
	for iRun in $(seq 1 $nRest); do
	
		# Bias correct mean EPI (less-aggressive N3).
		
		N3BiasFieldCorrection 3 \
			"${PWD}/meanEPI/meanRestEPI_run-$(printf %02d ${iRun}).nii" \
			"${PWD}/meanEPI/meanRestN3_run-$(printf %02d ${iRun}).nii.gz"
	
		meanRestEPI="${PWD}/meanEPI/meanRestN3_run-$(printf %02d ${iRun}).nii.gz"
		
		# Quick brain mask estimation.
		
		mri_synthstrip -i ${meanRestEPI} -m "${PWD}/meanEPI/meanRestN3_run-$(printf %02d ${iRun})_mask.nii.gz"
		ImageMath 3 "${PWD}/meanEPI/meanRestN3_run-$(printf %02d ${iRun})_mask_dil.nii.gz" MD "${PWD}/meanEPI/meanRestN3_run-$(printf %02d ${iRun})_mask.nii.gz" 10
		
		coregMaskEPI="${PWD}/meanEPI/meanRestN3_run-$(printf %02d ${iRun})_mask_dil.nii.gz"
		
		# Estimate coregistration.
		
		antsRegistration -d 3 \
			-o "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_" \
			-r [${ssT1},${meanRestEPI},0] \
			-t Translation[0.10] \
			-m Mattes[${ssT1},${meanRestEPI},1,32,Regular,0.30] \
			-c [12000x12000x11110,1e-8,20] \
			-s 4x2x1vox \
			-f 6x4x2 \
			-u 0 \
			-t Rigid[0.10] \
			-m Mattes[${ssT1},${meanRestEPI},1,32,Regular,0.30] \
			-c [12000x12000x11110,1e-8,20] \
			-s 4x2x1vox \
			-f 3x2x1 \
			-u 0 \
			-t Affine[0.10] \
			-m Mattes[${ssT1},${meanRestEPI},1,32,Regular,0.30] \
			-c [12000x12000x11110,1e-8,20] \
			-s 4x2x1vox \
			-f 3x2x1 \
			-u 0 \
			-t SyN[0.20,3,0] \
			-m Mattes[${ssT1},${meanRestEPI},0.5,32] \
			-m CC[${ssT1},${meanRestEPI},0.5,4] \
			-c [100x100x50,-0.01,5] \
			-s 1x0.5x0vox \
			-f 4x2x1 \
			-u 0 \
			-w [0.025,0.975] \
			-z 1 \
			-v 1
				
	done
	
fi

# Register T1 to MNI space.
#-----------------------------------------------------------------------------------------
# Specify smoothed segmentation probability maps for MNI template.

mniPost1="/zwork/tyler/splitMRI/template/MNI152_T1_2mm_spost1.nii.gz"
mniPost2="/zwork/tyler/splitMRI/template/MNI152_T1_2mm_spost2.nii.gz"
mniPost3="/zwork/tyler/splitMRI/template/MNI152_T1_2mm_spost3.nii.gz"
mniPost4="/zwork/tyler/splitMRI/template/MNI152_T1_2mm_spost4.nii.gz"
mniPost5="/zwork/tyler/splitMRI/template/MNI152_T1_2mm_spost5.nii.gz"
mniPost6="/zwork/tyler/splitMRI/template/MNI152_T1_2mm_spost6.nii.gz"

# And now for the patient.

subPost1="${pathToT1}/${patID}_T1w_smoothedPosteriors1.nii.gz"
subPost2="${pathToT1}/${patID}_T1w_smoothedPosteriors2.nii.gz"
subPost3="${pathToT1}/${patID}_T1w_smoothedPosteriors3.nii.gz"
subPost4="${pathToT1}/${patID}_T1w_smoothedPosteriors4.nii.gz"
subPost5="${pathToT1}/${patID}_T1w_smoothedPosteriors5.nii.gz"
subPost6="${pathToT1}/${patID}_T1w_smoothedPosteriors6.nii.gz"

# Grab the MNI template itself.

mniTemplate="/zwork/tyler/splitMRI/template/MNI152_T1_2mm_brain.nii.gz"

# Estimate registration with additional objectives at the nonlinear BSplineSyN stage.

antsRegistration -d 3 \
	-o "${pathToT1}/anat2mni_" \
	-r [${mniTemplate},${ssT1},1] \
	-t Translation[0.10] \
	-m Mattes[${mniTemplate},${ssT1},1,32,Regular,0.30] \
	-c [12000x12000x11110,1e-8,20] \
	-s 4x2x1vox \
	-f 6x4x2 \
	-u 0 \
	-t Rigid[0.10] \
	-m Mattes[${mniTemplate},${ssT1},1,32,Regular,0.30] \
	-c [12000x12000x11110,1e-8,20] \
	-s 4x2x1vox \
	-f 3x2x1 \
	-u 0 \
	-t Affine[0.10] \
	-m Mattes[${mniTemplate},${ssT1},1,32,Regular,0.30] \
	-c [12000x12000x11110,1e-8,20] \
	-s 4x2x1vox \
	-f 3x2x1 \
	-u 0 \
	-t BSplineSyN[0.15,26,0,3] \
	-m Mattes[${mniTemplate},${ssT1},1,32] \
	-m CC[${mniTemplate},${ssT1},1,4] \
	-m MSQ[${mniPost1},${subPost1},0.50,1] \
	-m MSQ[${mniPost2},${subPost2},0.75,1] \
	-m MSQ[${mniPost3},${subPost3},0.50,1] \
	-m MSQ[${mniPost4},${subPost4},0.75,1] \
	-m MSQ[${mniPost5},${subPost5},0.75,1] \
	-m MSQ[${mniPost6},${subPost6},0.75,1] \
	-c [100x100x50,-0.01,5] \
	-s 1x0.5x0vox \
	-f 4x2x1 \
	-x ${costFunMask} \
	-u 0 \
	-w [0.005,0.995] \
	-z 1 \
	-v 1

# Apply transforms.
#-----------------------------------------------------------------------------------------
# Warp T1 --> MNI space.

antsApplyTransforms -d 3 \
	-r ${mniTemplate} \
	-i ${ssT1} \
	-n LanczosWindowedSinc \
	-o "${pathToT1}/T1Warped.nii.gz" \
	-t "${pathToT1}/anat2mni_1Warp.nii.gz" \
	-t "${pathToT1}/anat2mni_0GenericAffine.mat" \
	-v 1

# Normalize whole-brain mask.

antsApplyTransforms -d 3 \
	-r ${mniTemplate} \
	-i ${brainMask} \
	-n GenericLabel \
	-o "${pathToT1}/BrainMaskWarped.nii.gz" \
	-t "${pathToT1}/anat2mni_1Warp.nii.gz" \
	-t "${pathToT1}/anat2mni_0GenericAffine.mat" \
	-v 1
	
# Erode WM and CSF masks by one voxel.

ImageMath 3 "${pathToT1}/csfEroded.nii.gz" ME "${pathToT1}/csfMask.nii.gz" 1
ImageMath 3 "${pathToT1}/wmEroded.nii.gz" ME "${pathToT1}/wmMask.nii.gz" 1	
	
# Finish additional transformations depending on the number of EPIs collected.

if [ $numEO -eq 1 ]; then
	
	# CSF mask + additional erosion.
	
	antsApplyTransforms -d 3 \
		-r ${meanRestEPI} \
		-i "${pathToT1}/csfEroded.nii.gz" \
		-n GenericLabel \
		-o "${PWD}/meanEPI/csfMaskRest.nii.gz" \
		-t ["${PWD}/meanEPI/mean2anat_rest_0GenericAffine.mat",1] \
		-t "${PWD}/meanEPI/mean2anat_rest_1InverseWarp.nii.gz" \
		-v 1
		
	ImageMath 3 "${PWD}/meanEPI/csfMaskRest.nii.gz" ME "${PWD}/meanEPI/csfMaskRest.nii.gz" 1
	
	# WM mask + additional erosion.
	
	antsApplyTransforms -d 3 \
		-r ${meanRestEPI} \
		-i "${pathToT1}/wmEroded.nii.gz" \
		-n GenericLabel \
		-o "${PWD}/meanEPI/wmMaskRest.nii.gz" \
		-t ["${PWD}/meanEPI/mean2anat_rest_0GenericAffine.mat",1] \
		-t "${PWD}/meanEPI/mean2anat_rest_1InverseWarp.nii.gz" \
		-v 1
		
	ImageMath 3 "${PWD}/meanEPI/wmMaskRest.nii.gz" ME "${PWD}/meanEPI/wmMaskRest.nii.gz" 1	

	# Obtain coregistered mean image (EPI --> T1).
		
	antsApplyTransforms -d 3 \
		-r ${ssT1} \
		-i ${meanRestEPI} \
		-n LanczosWindowedSinc \
		-o "${PWD}/meanEPI/meanRestCoreg.nii.gz" \
		-t "${PWD}/meanEPI/mean2anat_rest_1Warp.nii.gz" \
		-t "${PWD}/meanEPI/mean2anat_rest_0GenericAffine.mat" \
		-v 1
		
	# Obtain final warped mean EPI in MNI space.
		
	antsApplyTransforms -d 3 \
		-r ${mniTemplate} \
		-i ${meanRestEPI} \
		-n LanczosWindowedSinc \
		-o "${PWD}/meanEPI/meanRestWarped.nii.gz" \
		-t "${pathToT1}/anat2mni_1Warp.nii.gz" \
		-t "${pathToT1}/anat2mni_0GenericAffine.mat" \
		-t "${PWD}/meanEPI/mean2anat_rest_1Warp.nii.gz" \
		-t "${PWD}/meanEPI/mean2anat_rest_0GenericAffine.mat" \
		-v 1
		
	# Warp the 4D resting-state data.
	
	data4D="${pathToEPI}/au${patID}_task-rest_bold.nii"
					
	antsApplyTransforms -d 3 \
		-r ${mniTemplate} \
		-i ${data4D} \
		-e 3 \
		-n LanczosWindowedSinc \
		-o "${pathToEPI}/au${patID}_restWarped.nii.gz" \
		-t "${pathToT1}/anat2mni_1Warp.nii.gz" \
		-t "${pathToT1}/anat2mni_0GenericAffine.mat" \
		-t "${PWD}/meanEPI/mean2anat_rest_1Warp.nii.gz" \
		-t "${PWD}/meanEPI/mean2anat_rest_0GenericAffine.mat" \
		-v 1
		
else

	for jRun in $(seq 1 $nRest); do

		meanRestEPI="${PWD}/meanEPI/meanRestN3_run-$(printf %02d ${iRun}).nii.gz"
		
		# CSF mask + additional erosion.
	
		antsApplyTransforms -d 3 \
			-r ${meanRestEPI} \
			-i "${pathToT1}/csfEroded.nii.gz" \
			-n GenericLabel \
			-o "${PWD}/meanEPI/csfMaskRest_run-$(printf %02d ${iRun}).nii.gz" \
			-t ["${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_0GenericAffine.mat",1] \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_1InverseWarp.nii.gz" \
			-v 1
			
		ImageMath 3 "${PWD}/meanEPI/csfMaskRest_run-$(printf %02d ${iRun}).nii.gz" ME "${PWD}/meanEPI/csfMaskRest_run-$(printf %02d ${iRun}).nii.gz" 1
		
		# WM mask + additional erosion.
		
		antsApplyTransforms -d 3 \
			-r ${meanRestEPI} \
			-i "${pathToT1}/wmEroded.nii.gz" \
			-n GenericLabel \
			-o "${PWD}/meanEPI/wmMaskRest_run-$(printf %02d ${iRun}).nii.gz" \
			-t ["${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_0GenericAffine.mat",1] \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_1InverseWarp.nii.gz" \
			-v 1
			
		ImageMath 3 "${PWD}/meanEPI/wmMaskRest_run-$(printf %02d ${iRun}).nii.gz" ME "${PWD}/meanEPI/wmMaskRest_run-$(printf %02d ${iRun}).nii.gz" 1	
	
		# Obtain coregistered mean image (EPI --> T1).
			
		antsApplyTransforms -d 3 \
			-r ${ssT1} \
			-i ${meanRestEPI} \
			-n LanczosWindowedSinc \
			-o "${PWD}/meanEPI/meanRestCoreg_run-$(printf %02d ${iRun}).nii.gz" \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_1Warp.nii.gz" \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_0GenericAffine.mat" \
			-v 1
			
		# Obtain final warped mean EPI in MNI space.
			
		antsApplyTransforms -d 3 \
			-r ${mniTemplate} \
			-i ${meanRestEPI} \
			-n LanczosWindowedSinc \
			-o "${PWD}/meanEPI/meanRestWarped_run-$(printf %02d ${iRun}).nii.gz" \
			-t "${pathToT1}/anat2mni_1Warp.nii.gz" \
			-t "${pathToT1}/anat2mni_0GenericAffine.mat" \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_1Warp.nii.gz" \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_0GenericAffine.mat" \
			-v 1
			
		# Warp the 4D resting-state data.
		
		data4D="${pathToEPI}/au${patID}_task-rest_run-$(printf %02d ${iRun})_bold.nii"
						
		antsApplyTransforms -d 3 \
			-r ${mniTemplate} \
			-i ${data4D} \
			-e 3 \
			-n LanczosWindowedSinc \
			-o "${pathToEPI}/au${patID}_restWarped_run-$(printf %02d ${iRun}).nii.gz" \
			-t "${pathToT1}/anat2mni_1Warp.nii.gz" \
			-t "${pathToT1}/anat/anat2mni_0GenericAffine.mat" \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_1Warp.nii.gz" \
			-t "${PWD}/meanEPI/mean2anat_rest_run-$(printf %02d ${iRun})_0GenericAffine.mat" \
			-v 1
			
	done
	
fi

exit 0