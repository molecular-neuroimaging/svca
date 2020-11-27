function [normGrayTAC, normWhiteTAC, normIDIF, normHighBindTAC] = superX_extractNormTACs(dynPET,brainMask, whiteMask, grayMask, highBindMask, time)
% Tissue classes include white matter, gray matter, blood, high binding gray matter
%
%   INPUTS
%         dynPET: 4D matrix containing dynamic PET data
%         brainMask: 3D matrix containing binary brain data specifying voxels of brain mask. Must be the same mask used for extracting the kinetic classes!
%         highBindMask: 3D matrix containing binary data specifying voxels of high binding tissue (ex. thalamus)
%         time: 1D double with PET mid-frame times (in minutes)
%   OUTPUTS
%         normTACs: (normGrayTAC, normWhiteTAC, normIDIF, normHighBindTAC) 2D double with normalised time activity curves in the
%         following order: GM, WM,image-derived blood input function, high binding gray matter

%Copyright - Julia Schubert and Mattia Veronese - London 2020
%-------------------------------------------------------------------------
%% GENERAL SETTINGS
normGrayTAC = [];
normWhiteTAC = [];
normIDIF=[];
normHighBindTAC=[];

%% LOAD DATA
%Reading PET
dynPET = double(dynPET);
sizePET = size(dynPET);
nFrame = sizePET(4);

if length(time) ~= nFrame %check consistency: time vector should match with the number of PET frames
    disp('ERROR: image dimensions of dynPET and time must agree');
    return;
end

%Creating variable containing the frame durations
delta=zeros(nFrame,1);
delta(1)=2*time(1);
for k=2:nFrame
    delta(k)=abs(2*(time(k)-sum(delta)));
end


%% PET DATA NORMALISATION - THIS IS VERY IMPORTANT brainMask needs to match the same one used for the extraction of the supervised reference region
fprintf('\nNormalising dynamic PET....\n');

[normDYN]=dynPET_normalization(dynPET, brainMask); %Normalise dynamic PET data using brainMask
%This is important - brainMask must be the same used for extracting
%the supervised reference region with the generated kinetic
%classes!

disp('Normalisation complete');


%% Extract normalised TACs (GM, WM and High Binding)
fprintf('\nExtracting normalised TACs....\n');

v               = zeros(nFrame,1);
normGrayTAC     = zeros(nFrame,1);
normWhiteTAC	= zeros(nFrame,1);
normHighBindTAC	= zeros(nFrame,1);

for z=1:sizePET(3)
    for x=1:sizePET(1)
        for y=1:sizePET(2)
            if (whiteMask(x,y,z)>0)
                v(:) 		= normDYN(x,y,z,:); %using non-normalized PET for TAC extraction
                normWhiteTAC		=  normWhiteTAC+v;
            end
            if (grayMask(x,y,z)>0)
                v(:) 		= normDYN(x,y,z,:);
                normGrayTAC		=  normGrayTAC+v;
            end
            if (highBindMask(x,y,z)>0)
                v(:) 		= normDYN(x,y,z,:);
                normHighBindTAC		=  normHighBindTAC+v; 
            end
            
        end
        
    end
end

ref_idx_white    = find(whiteMask>0);
ref_idx_gray     = find(grayMask>0);
ref_idx_highBind = find(highBindMask>0);

normWhiteTAC    = normWhiteTAC/length(ref_idx_white);
normGrayTAC		=  normGrayTAC/length(ref_idx_gray);
normHighBindTAC = normHighBindTAC/length(ref_idx_highBind);

%% Extraction of normalised blood TAC
timeToPeak =1.5; % Time window within which blood peak is expected (in mins)
normIDIF = computeIDIF(normDYN,transpose(delta),timeToPeak, brainMask);

fprintf('Normalised TAC extraction complete\n\n\n');

end


