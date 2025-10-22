function splits_corePreproc(patID)
% Preprocess an individual patient's data.
%
% FORMAT splits_corePreproc(subID)
%
% REQUIRED INPUT:
%   patID
%       String/character array specifying patient ID.
%__________________________________________________________________________
%
% This wrapper performs the following preprocessing steps:
%        1) Head motion correction.
%        2) Slicetime correction.
%        3) Anatomical preprocessing.
%        4) Spatial registration.
%        5) Spatial smoothing.
%        6) Anatomical/physiological noise estimation.
%        7) Global intensity scaling.
%        8) Linear detrending.
%        9) Nuisance signal regression.
%       10) Regional signal extraction.
%       11) Bandpass filtering.
%__________________________________________________________________________
%
% Part of the preprocessing and analysis code for:
% Santander et al. (2025). Full interhemispheric integration sustained by a
% fraction of posterior callosal fibers. Proc Natl Acad Sci USA, 122(43),
% e2520190122. https://doi.org/10.1073/pnas.2520190122
%
% Github repo: https://github.com/tsantander/splitBrainNetworks
%__________________________________________________________________________

% Preliminary setup.
%--------------------------------------------------------------------------
addpath('/zwork/tyler/splitMRI/code');
addpath('/home/tyler/Documents/MATLAB/spm/spm12');
addpath(genpath('/home/tyler/Documents/MATLAB/BrainGraphs'));
system('export ANTSPATH=/sw/ANTs/bin');
system('export PATH=${ANTSPATH}:${PATH}');

% Specify the parent directory for the patient data, navigate to it.

parentDir = '/zwork/tyler/splitMRI/data';
cd(parentDir);

% Run initial processing pipeline.
%--------------------------------------------------------------------------
% Start the clock so we can track computation time.

procStart = tic;
disp(' ');
disp(['|| Running ' patID '. Please wait...']);

% Jump into this patient's data directory.

cd([pwd '/' patID]);

% Get number of resting-state scans collected for this patient.

switch patID
    case {'sub-01', 'sub-04', 'sub-05'}
        nRest = 1;
    case {'sub-02', 'sub-06'}
        nRest = 2;
    case 'sub-07'
        nRest = 3;
end

% Motion correction and slicetime correction.

splits_realignUnwarp;
splits_sliceTimeCorrect;

% Anatomical preprocessing.

system(['/zwork/tyler/splitMRI/code/splits_processAnat.sh ' patID]);

% Coregistration and spatial normalization using cost-function masking.

switch patID
    case 'sub-01'
        cfMask = '/zwork/tyler/splitMRI/template/cfmAntSplit.nii.gz';
    case 'sub-02'
        cfMask = '/zwork/tyler/splitMRI/template/cfmSplOnly.nii.gz';
    case {'sub-04', 'sub-05', 'sub-06', 'sub-07'}
        cfMask = '/zwork/tyler/splitMRI/template/cfmFullSplit.nii.gz';
end

system(['/zwork/tyler/splitMRI/code/splits_spatialRegCFM.sh ' patID ' ' cfMask ' ' num2str(nRest)]);

% Smooth data and estimate anatomical noise components.

splits_smooth;
splits_estimateAnatomicalNoise(nRest);

% Run additional denoising necessary for rsFC analyses.
%--------------------------------------------------------------------------
for iRest = 1:nRest

    if (nRest == 1)
        basename = ['sau' patID '_restWarped.nii.gz'];
    else
        basename = ['sau' patID '_restWarped_' sprintf('run-%02d', iRest) '.nii.gz'];
    end

    % Scale data to a global median = 1000.

    brainMask = [pwd '/anat/BrainMaskWarped.nii.gz'];

    bgt_globalNorm([pwd '/func/' basename], 'median', brainMask);

    % Estimate and remove linear trend.

    bgt_detrend([pwd '/func/g' basename], 1);

    if (nRest == 1)
        unix(['mv -f ' pwd '/func/dtMatrix.mat ' pwd '/func/dtMatrixRest.mat']);
    else
        unix(['mv -f ' pwd '/func/dtMatrix.mat ' pwd '/func/dtMatrixRest_' sprintf('run-%02d', iRest) '.mat']);
    end

    % Detrend motion parameters / anatomical noise components.

    if (nRest == 1)
        
        motionParams = dlmread([pwd '/func/' spm_select('List', [pwd '/func'], '^rp.*rest.*txt$')]);

        load([pwd '/func/wmNoiseRest.mat'], 'wmNoise');
        load([pwd '/func/csfNoiseRest.mat'], 'csfNoise');
        load([pwd '/func/globalSignalRest.mat'], 'globalSignal');
        load([pwd '/func/dtMatrixRest.mat'], 'R');

    else

        motionParams = dlmread([pwd '/func/' spm_select('List', [pwd '/func'], ['^rp.*rest_' sprintf('run-%02d', iRest) '.*txt$'])]);

        load([pwd '/func/wmNoiseRest_' sprintf('run-%02d', iRest) '.mat'], 'wmNoise');
        load([pwd '/func/csfNoiseRest_' sprintf('run-%02d', iRest) '.mat'], 'csfNoise');
        load([pwd '/func/globalSignalRest_' sprintf('run-%02d', iRest) '.mat'], 'globalSignal');
        load([pwd '/func/dtMatrixRest_' sprintf('run-%02d', iRest) '.mat'], 'R');

    end

    dtMotion    = R*motionParams;
    dtWM        = R*wmNoise;
    dtCSF       = R*csfNoise;
    dtGlobalSig = R*globalSignal;

    % Concatenate detrended nuisance regressors for 36P model.

    allNuisance = [dtMotion, dtWM(:,1), dtCSF(:,1), dtGlobalSig];

    % Check to see if we need to model additional spike regressors.

    dtMotion(:,4:6) = dtMotion(:,4:6) .* 50;
    ddtMP           = diff(dtMotion);
    fwd             = [0; sum(abs(ddtMP), 2)];

    spikeID = find(fwd >= 0.50);

    if (~isempty(spikeID))

        spikeReg = zeros(size(fwd,1), length(spikeID));

        for iSpike = 1:length(spikeID)

            spikeReg(spikeID(iSpike), iSpike) = 1;

        end

    else

        spikeReg = [];

    end 

    % Do confound regression (spike regressors will be ignored if empty).

    bgt_regressNuisance([pwd '/func/dg' basename], allNuisance, 'expandDeriv', spikeReg);

    % Extract regional timeseries from 432 node composite (Schaefer/Tian).

    atlasImg = '/zwork/tyler/splitMRI/template/compositeAtlas_7Network_432Node.nii.gz';

    [timeSeries] = bgt_extractRegionalTimeseries([pwd '/func/ndg' basename], atlasImg, 'eigen1');

    if (nRest == 1)
        unix(['mv -f ' pwd '/func/timeSeries.mat ' pwd '/func/timeSeriesRest.mat']);
    else
        unix(['mv -f ' pwd '/func/timeSeries.mat ' pwd '/func/timeSeriesRest_' sprintf('run-%02d', iRest) '.mat']);
    end

    % Bandpass filter regional timeseries.

    bgt_bandpass(timeSeries, 2, [0.01, 0.08], 'butter');
    
    if (nRest == 1)
        unix(['mv -f ' pwd '/func/bandpassSeries.mat ' pwd '/func/bandpassSeriesRest.mat']);
    else
        unix(['mv -f ' pwd '/func/bandpassSeries.mat ' pwd '/func/bandpassSeriesRest_' sprintf('run-%02d', iRest) '.mat']);
    end

end

% All done!
%--------------------------------------------------------------------------
% Display computation time.

procEnd = toc(procStart);
disp(['|| Done with core preprocessing for ' patID ' in ' num2str(procEnd/60) ' minutes']);
clear procStart procEnd

% Return to parent directory.

cd(parentDir);
        
end
