%% Schubert Julia - Test Script for superX_extractNormTACs.m

% Initialization
clc;
close all;

%Change current directory to code directory
if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

%Add utility folders to path
codeDir = pwd;
dataDir = fullfile(fileparts(codeDir), 'data');
addpath(genpath(codeDir),'-begin');

%% MAIN
%Initialize directory cycle for control subjects (i.e TRAINING DATASET)
cd(dataDir);
listSubs=dir('exampleCON*');

for s=1:size(listSubs,1)
    %% Define current subject
    subjID=listSubs(s).name;
    cd(subjID);
    fprintf(['Processing subject: ',num2str(subjID), '\n']);
    
    %% LOAD INPUT
    %Read PET file
    dynPET_file = [subjID,'_dynPETforSuperRef.nii'];
    if (exist([dynPET_file, '.gz'], 'file'))
        gunzip([dynPET_file, '.gz']);
    end
    if (~exist([subjID,'_dynPETforSuperRef.nii'], 'file'))
        fprintf(['File missing: ', dynPET_file, '\n']);
        continue;
    end
    
    niiPET=load_nii(dynPET_file);
    dynPET=double(niiPET.img)/1000; %from Bq to KBq
    sizePET=size(dynPET);
    
    %Load mask for gray matter tissue class (use only cortex and cerebellum)
    GM_file = ([subjID, '_grayMapForSuperRef.nii']);
    niiGM = load_nii(GM_file);
    GM = double(niiGM.img);
    GM = GM./max(GM(:)); %Normalize GM between 0 and 1
    
    innerGM_file = ([subjID, '_innerGMForSuperRef.nii']);
    niiInnerGM = load_nii(innerGM_file);
    innerGM = double(niiInnerGM.img);
    indInnerGM = find(innerGM>0.9);
    
    indMaskExtraction=find(GM>0.9);
    grayMask = zeros(sizePET(1:3));
    grayMask(indMaskExtraction) = 1;
    grayMask(indInnerGM) = 0;
    
    %Load mask for white matter tissue class
    WM_file = ([subjID, '_whiteMapForSuperRef.nii']);
    niiWM = load_nii(WM_file);
    WM = double(niiWM.img);
    WM = WM./max(WM(:)); %Normalize WM between 0 and 1
    
    indMaskWhite=find(WM>0.9);
    refMaskWhite = zeros(sizePET(1:3));
    whiteMask = zeros(sizePET(1:3));
    whiteMask(indMaskWhite) = 1;
    
    %Load mask for high binding tissue class
    thalamus_file = ([subjID, '_thalamusMaskForSuperRef.nii']);
    niiThalamus = load_nii(thalamus_file);
    highBindMask = double(niiThalamus.img);
    
    
    %Create whole brain mask for normalization of dynamic PET
    BRAIN_file = ([subjID, '_brainMaskForSuperRef.nii']);
    brainMask = load_nii(BRAIN_file);
    brainMask = double(brainMask.img);
    
    %Load PET mid-frame time
    load('time'); %PLEASE MAKE SURE IT IS IN MINUTES
    
    %RUN superX_extractNormTACs
    [normGrayTAC, normWhiteTAC, normIDIF, normHighBindTAC] = superX_extractNormTACs(dynPET,brainMask, whiteMask, grayMask, highBindMask, time);
    
    %PLOT INDIVIDUAL NORMALISED TACs (PLEASE SAVE THEM)
    figure
    plot(time,[normGrayTAC, normWhiteTAC, normIDIF, normHighBindTAC]);
    xlabel('TIME (mins)')
    ylabel('NORMSALISED ACTIVITY')
    legend('GM','WM','Blood','HighBindGM')
    title(subjID)
    
    %EXIT the folder
    cd(dataDir);
    
end

% Remove added paths
rmpath(genpath(codeDir));



