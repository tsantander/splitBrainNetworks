function splits_sliceTimeCorrect
% Slice-time correction.
%
% FORMAT splits_sliceTimeCorrect
%__________________________________________________________________________
%
% Perform slice-time correction (assumes this script is being called from a 
% wrapper and we are already in a patient's base data directory).
%__________________________________________________________________________
%
% Part of the preprocessing and analysis code for:
% Santander et al. (2025). Full interhemispheric integration sustained by a
% fraction of posterior callosal fibers. Proc Natl Acad Sci USA, 122(43),
% e2520190122. https://doi.org/10.1073/pnas.2520190122
%
% Github repo: https://github.com/tsantander/splitBrainNetworks
%__________________________________________________________________________

% Slice-time correct resting-state data.
%--------------------------------------------------------------------------
% Find all resting-state scans for this patient (following motion
% correction) and get corresponding metadata files as well.

restEPI  = cellstr(strcat([pwd '/func/'], spm_select('List', [pwd '/func'], '^usub.*task-rest.*bold.nii$')));
restJSON = cellstr(strcat([pwd '/func/'], spm_select('List', [pwd '/func'], '^sub.*task-rest.*bold.json$')));

% Now loop over scans (if necessary).

for iRest = 1:length(restEPI)
    
    setDefaultsSPM;

    matlabbatch = {};

    % Specify scans.

    matlabbatch{1}.spm.temporal.st.scans = {restEPI(iRest)};

    % Read in json.

    fid     = fopen(restJSON{iRest});
    raw     = fread(fid,inf);
    str     = char(raw');
    metaDat = jsondecode(str);
    fclose(fid);

    % Specify STC parameters and run.

    matlabbatch{1}.spm.temporal.st.nslices  = length(metaDat.SliceTiming);
    matlabbatch{1}.spm.temporal.st.tr       = 2;
    matlabbatch{1}.spm.temporal.st.ta       = 0;
    matlabbatch{1}.spm.temporal.st.so       = metaDat.SliceTiming * 1000;
    matlabbatch{1}.spm.temporal.st.refslice = 0;
    matlabbatch{1}.spm.temporal.st.prefix   = 'a';

    spm_jobman('run', matlabbatch);
    
end

end

%-------------------------------------------------------------------------%
% BEGIN SUBROUTINES                                                       %
%-------------------------------------------------------------------------%

% Initialize default parameters for SPM.
%-------------------------------------------------------------------------%
function setDefaultsSPM

spm('defaults','fMRI');
spm_jobman('initcfg');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% END SUBROUTINES                                                         %
%-------------------------------------------------------------------------%