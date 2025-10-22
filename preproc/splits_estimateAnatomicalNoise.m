function splits_estimateAnatomicalNoise(nRest)
% Estimate anatomical/physiological noise components.
%
% FORMAT splits_estimateAnatomicalNoise
%__________________________________________________________________________
%
% This function will estimate components of anatomical/physiological noise
% for later nuisance signal regression (via SVD). We'll also obtain mean
% global signal at each TR. This script assumes it is being called from a 
% wrapper and we have already changed into a patient's base data directory).
%__________________________________________________________________________
%
% Part of the preprocessing and analysis code for:
% Santander et al. (2025). Full interhemispheric integration sustained by a
% fraction of posterior callosal fibers. Proc Natl Acad Sci USA, 122(43),
% e2520190122. https://doi.org/10.1073/pnas.2520190122
%
% Github repo: https://github.com/tsantander/splitBrainNetworks
%__________________________________________________________________________

% Get resting-state noise estimates.
%--------------------------------------------------------------------------
% Specify where the relevant data live.

anatPath = [pwd '/anat'];
meanPath = [pwd '/meanEPI']; 
funcPath = [pwd '/func'];

% Loop over runs (if applicable).
        
for iRun = 1:nRest

    % Load in HMC + STC data (prior to spatial normalization).

    if (nRest == 1)
        data4D  = load_nii(strcat([funcPath '/'], spm_select('List', funcPath, '^ausub.*task-rest.*bold.nii$')));
    else
        data4D  = load_nii(strcat([funcPath '/'], spm_select('List', funcPath, ['^ausub.*task-rest_' sprintf('run-%02d', iRun) '.*bold.nii$'])));
    end
        
    [x,y,z,t] = size(data4D.img);
    data2D    = reshape(data4D.img, x*y*z, t);
    data2D    = double(data2D)';

    clear data4D

    % Load in WM mask in same EPI space - get noise components.
    
    if (nRest == 1)
        mask3D = load_nii([meanPath '/wmMaskRest.nii.gz']);
    else
        mask3D = load_nii([meanPath '/wmMaskRest_' sprintf('run-%02d', iRun) '.nii.gz']);
    end

    mask2D = reshape(mask3D.img, 1, numel(mask3D.img));
    
    clear mask3D

    disp('|| Extracting first 5 principal components from WM mask');

    noiseData = data2D(:, logical(mask2D));
    zNoise    = zscore(noiseData); 
    wmNoise   = projectComponents(zNoise);

    disp('|| Finished component extraction. Saving...');

    if (nRest == 1)
        save([funcPath '/wmNoiseRest.mat'], 'wmNoise');
    else
        save([funcPath '/wmNoiseRest_' sprintf('run-%02d', iRun) '.mat'], 'wmNoise');
    end

    clear noiseData zNoise wmNoise

    % Load CSF mask in EPI space - get noise components.

    if (nRest == 1)
        mask3D = load_nii([meanPath '/csfMaskRest.nii.gz']);
    else
        mask3D = load_nii([meanPath '/csfMaskRest_' sprintf('run-%02d', iRun) '.nii.gz']);
    end

    mask2D = reshape(mask3D.img, 1, numel(mask3D.img));

    clear mask3D

    disp('|| Extracting first 5 principal components from CSF mask');

    noiseData = data2D(:, logical(mask2D));
    zNoise    = zscore(noiseData); 
    csfNoise  = projectComponents(zNoise);

    disp('|| Finished component extraction. Saving...');

    if (nRest == 1)
        save([funcPath '/csfNoiseRest.mat'], 'csfNoise');
    else
        save([funcPath '/csfNoiseRest_' sprintf('run-%02d', iRun) '.mat'], 'csfNoise');
    end

    clear noiseData zNoise csfNoise data2D

    % Now load in smoothed/normalized resting-state data.

    if (nRest == 1)
        data4D = load_nii(strcat([funcPath '/'], spm_select('List', funcPath, '^sausub.*restWarped.nii.gz$')));
    else
        data4D = load_nii(strcat([funcPath '/'], spm_select('List', funcPath, ['^sausub.*restWarped_' sprintf('run-%02d', iRun) '.nii.gz$'])));
    end

    % Get brain mask in MNI space and toss non-brain voxels.

    mask3D = load_nii([anatPath '/BrainMaskWarped.nii.gz']);
    mask2D = reshape(mask3D.img, 1, numel(mask3D.img));

    [x,y,z,t]                   = size(data4D.img);
    data2D                      = reshape(data4D.img, x*y*z, t);
    data2D(~logical(mask2D),:)  = [];

    clear data4D mask3D mask2D

    % Get global signal vector.

    disp('|| Estimating global signal over time');

    globalSignal = mean(data2D)';

    disp('|| Finished. Saving...');

    if (nRest == 1)
        save([funcPath '/globalSignalRest.mat'], 'globalSignal');
    else
        save([funcPath '/globalSignalRest_' sprintf('run-%02d', iRun) '.mat'], 'globalSignal');
    end

    clear data2D mask2D globalSignal

end

end

%-------------------------------------------------------------------------%
% BEGIN SUBROUTINES                                                       %
%-------------------------------------------------------------------------%

% Project noise data into component space.
%--------------------------------------------------------------------------
function [anatNoise] = projectComponents(zNoise)

[m,n]     = size(zNoise);
anatNoise = zeros(m,5);

if (m > n)
    
    [~, eVal, eVect] = svd(zNoise'*zNoise);
    eVal             = diag(eVal);
    
    for iComp = 1:5
        
        v                  = eVect(:,iComp);
        u                  = zNoise*v/sqrt(eVal(iComp));

        direction          = sign(sum(v));
        b                  = u*direction;
        anatNoise(:,iComp) = b*sqrt(eVal(iComp)/n);
        
    end
    
else
    
    [~, eVal, eVect] = svd(zNoise*zNoise');
    eVal             = diag(eVal);
    
    for iComp = 1:5
        
        u                  = eVect(:,iComp);
        v                  = zNoise'*u/sqrt(eVal(iComp));

        direction          = sign(sum(v));
        b                  = u*direction;
        anatNoise(:,iComp) = b*sqrt(eVal(iComp)/n);
        
    end

end

end
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------%
% END SUBROUTINES                                                         %
%-------------------------------------------------------------------------%