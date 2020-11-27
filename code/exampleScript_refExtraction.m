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
%Initialize directory cycle for patients (i.e TESTING DATASET)
cd(dataDir);
listSubs=dir('examplePAT*');

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
    
    DIM=size(dynPET);   %PET dimension
    lt=DIM(4);    %number of frames
    
    %Read kinetic classes
    file_kineticClasses = fullfile(superRefUtility, 'tissueClassesForSupervisedRef.csv');
    
    %Load mask for gray matter tissue class (using only cortex and cerebellum)
    GM_file = ([subjID, '_grayMapForSuperRef.nii']);
    niiGM = load_nii(GM_file);
    GM = double(niiGM.img);
    GM = GM./max(GM(:)); %Normalize GM between 0 and 1
    
    indMaskExtraction=find(GM>0.9);
    refMaskExtraction = zeros(DIM(1:3));
    refMaskExtraction(indMaskExtraction) = 1;
    
    %Create whole brain mask for normalization of dynamic PET
    BRAIN_file = ([subjID, '_brainMaskForSuperRef.nii']);
    brainMask = load_nii(BRAIN_file);
    brainMask = double(brainMask.img);
    
    %Load PET mid-frame time
    load('time'); %PLEASE MAKE SURE IT IS IN MINUTES
    
    %RUN superX_referExtraction
    [RefTAC, RefMask] = superX_referExtraction(dynPET,file_kineticClasses,brainMask,refMaskExtraction, time);
    
    %PLOT INDIVIDUAL NORMALISED TACs (PLEASE SAVE THEM)
    if ~isempty(RefTAC)
        figure
        plot(time,RefTAC);
        xlabel('TIME (mins)')
        ylabel('Supervised reference region')
        title(subjID)
    end
    
    %EXIT the folder
    cd(dataDir);
end

% Remove added paths
rmpath(genpath(codeDir));



