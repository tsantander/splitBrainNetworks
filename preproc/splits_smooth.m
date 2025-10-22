function splits_smooth
% Spatial smoothing.
%
% FORMAT splits_smooth
%__________________________________________________________________________
%
% Applies a 5mm^3 FWHM smoothing kernel to the 4D restng-state data 
% (assumes this script is being called from a wrapper and we have already
% changed into a patient's base data directory).
%__________________________________________________________________________
%
% Part of the preprocessing and analysis code for:
% Santander et al. (2025). Full interhemispheric integration sustained by a
% fraction of posterior callosal fibers. Proc Natl Acad Sci USA, 122(43),
% e2520190122. https://doi.org/10.1073/pnas.2520190122
%
% Github repo: https://github.com/tsantander/splitBrainNetworks
%__________________________________________________________________________

% Initialize default SPM configurations for fMRI.
%--------------------------------------------------------------------------
setDefaultsSPM;

% Get all the normalized, realigned/unwarped functional data.
%--------------------------------------------------------------------------
% Might be multiple resting-state scans so check for that and gunzip.

rest4D = cellstr(strcat([pwd '/func/'], spm_select('List', [pwd '/func'], '^ausub.*restWarped.*nii.gz$')));

for iRest = 1:length(rest4D)
    unix(['gunzip -f ' rest4D{iRest}]);
end

% Now grab all the gunzipped files.

rest4D = cellstr(strcat([pwd '/func/'], spm_select('List', [pwd '/func'], '^ausub.*restWarped.*nii$')));

% Specify smoothing parameters.
%--------------------------------------------------------------------------
matlabbatch = {};

matlabbatch{1}.spm.spatial.smooth.data   = rest4D;
matlabbatch{1}.spm.spatial.smooth.fwhm   = [5 5 5];
matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
matlabbatch{1}.spm.spatial.smooth.im     = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

% Run spatial smoothing.
%--------------------------------------------------------------------------
spm_jobman('run',matlabbatch);

% Compress outputs.
%--------------------------------------------------------------------------
rest4D = cellstr(strcat([pwd '/func/'], spm_select('List', [pwd '/func'], '^ausub.*restWarped.*nii$')));

for iRest = 1:length(rest4D)
    unix(['gzip -f ' rest4D{iRest}]);
end

rest4D = cellstr(strcat([pwd '/func/'], spm_select('List', [pwd '/func'], '^sausub.*restWarped.*nii$')));

for iRest = 1:length(rest4D)
    unix(['gzip -f ' rest4D{iRest}]);
end

end

%-------------------------------------------------------------------------%
% BEGIN SUBROUTINES                                                       %
%-------------------------------------------------------------------------%

% Initialize default parameters for SPM.
%--------------------------------------------------------------------------
function setDefaultsSPM

spm('defaults','fMRI');
spm_jobman('initcfg');

end
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------%
% END SUBROUTINES                                                         %
%-------------------------------------------------------------------------%