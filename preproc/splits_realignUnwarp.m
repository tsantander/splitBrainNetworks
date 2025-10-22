function splits_realignUnwarp
% Motion realignment and unwarping.
%
% FORMAT splits_realignUnwarp
%__________________________________________________________________________
%
% Perform head motion correction and motion x susceptibility correction
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

% Realign/unwarp resting-state.
%--------------------------------------------------------------------------
% Find all resting-state scans for this patient.

restEPIs = cellstr(strcat([pwd '/func/'], spm_select('List', [pwd '/func'], '^sub.*task-rest.*bold.nii$')));

% Loop over if necessary.

for iRest = 1:length(restEPIs)
    
    % Initialize.

    [matlabbatch] = setParameters;

    % Grab the data.

    matlabbatch{1}.spm.spatial.realignunwarp.data.scans  = cellstr(restEPIs{iRest});
    matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = {};

    % Run the job.

    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
end

% Create new folder for the mean functional image(s) and move them there.
%--------------------------------------------------------------------------
mkdir([pwd '/meanEPI']);

% Check if this patient has multiple resting-state runs and handle
% accordingly.

if (length(restEPIs) > 1)

    for iRest = 1:length(restEPIs)

        [~,name,~] = fileparts(restEPIs{iRest});

        unix(['mv ' pwd '/func/meanu' name '.nii ' pwd '/meanEPI/meanRestEPI_' sprintf('run-%02d', iRest) '.nii']);

    end
    
else
    
    unix(['mv ' pwd '/func/meanu*rest*.nii ' pwd '/meanEPI/meanRestEPI.nii']);
    
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

% Define all parameters for realign/unwarp estimation.
%--------------------------------------------------------------------------
function [matlabbatch] = setParameters

matlabbatch = {};

% Set parameters for realignment stage.

matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 1;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep     = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm    = 5;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm     = 0;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 7;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap   = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight  = '';

% Set parameters for unwarpgin stage.

matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn   = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda   = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm       = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot      = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot      = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm   = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem      = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi      = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';

% Set parameters for reslicing.

matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap    = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask    = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix  = 'u';

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% END SUBROUTINES                                                         %
%-------------------------------------------------------------------------%