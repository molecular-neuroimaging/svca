function [RefTAC,RefMask ] = superX_referExtraction(dynPET,kineticClasses,brainMask,refMaskExtraction, time)
%superX_referenceExtraction extracts the reference region for a tracer X using
%supervised clustering given a set of predefined kinetic classes 
%   INPUTS
%         dynPET: 4D matrix containing dynamic PET data
%         file_kineticClasses: path to .csv file containing kinetic classes. Columns in .csv must include in this order: time, gray matter, white matter, blood, high binding gray matter
%         brainMask: 3D matrix containing binary brain data specifying voxels of brain mask. Must be the same mask used for extracting the kinetic classes!
%         refMaskExtraction: 3D matrix containing binary data specifying voxels to select from for supervised reference region
%         time: 1D double with PET mid-frame times
%   OUTPUTS
%         RefTAC: vector containing time activity curve of extracted supervised reference region
%         RefMask: 3D matrix showing voxels extracted for supervised reference region

%Copyright - Julia Schubert and Mattia Veronese - London 2020
%------------------------------------------------------------

%% GENERAL SETTINGS
threshold = 0.9; %Probability threshold to be reference region - based on literature
RefTAC =[];
RefMask =[];

%% LOAD DATA
%reading kinetic classes
a = builtin('load',kineticClasses);

[m,n]	= size(a);
tref	= a(:,1);
class	= a(:,2:n);

%reading PET
dynPET = double(dynPET);
sizePET = size(dynPET);
nFrame = sizePET(4);

%reading MASKs
sizeBRAIN = size(brainMask);
sizeREF = size(refMaskExtraction);

if sizePET(1:3) ~= sizeBRAIN
    disp('ERROR: image dimensions of dynPET and brainMask must agree');
    return;
end
if sizePET(1:3) ~= sizeREF
    disp('ERROR: image dimensions of dynPET and refMaskExtraction must agree');
    return;
end
if length(time) ~= nFrame %check consistency: time vector should match with the number of PET frames
    disp('ERROR: image dimensions of dynPET and time must agree');
    return;
end

%% PET DATA NORMALISATION - THIS IS VERY IMPORTANT brainMask needs to match the same one used for the generation of kinetic classes
fprintf('\nNormalising dynamic PET....\n');
[normDYN]=dynPET_normalization(dynPET, brainMask);
disp('Normalisation complete');

%% TIME INTERPOLATION
fprintf('\nInterpolating reference classes on dynamic PET....\n');
for k=nFrame:-1:1 %Select the last time frame from dynPET to be used for extraction of reference region
    if time(k)<=tref(m)
        break;
    end
end

maxk	= k;
tref	= [0 tref']';
class	= [zeros(n-1,1) class']';
classI	= zeros(maxk,n-1);

%interpolation of kinetic classes from original time to dynPET time
for k=1:n-1
    classI(:,k) = interp1(tref,class(:,k),time(1:maxk));
end
disp('Interpolation complete');

%% REFERENCE EXTRACTION
fprintf('\nExtracting supervised reference region....');
GRAY  = zeros(sizePET(1:3));
GRAYRATIO  = zeros(sizePET(1:3));
v		= zeros(nFrame,1);
frac=0.1;

for z=1:sizePET(3)
    for x=1:sizePET(1)
        for y=1:sizePET(2)
            if(brainMask(x,y,z)>0)
                v(:) 		= normDYN(x,y,z,1:maxk); %using normalized PET for class identification
                bv          = lsqnonneg(classI,v);
                GRAY(x,y,z) 	= bv(1);
                if (sum(bv)>0)
                    GRAYRATIO(x,y,z) = GRAY(x,y,z)/sum(bv);
                end
            end
        end
    end
    if (z/sizePET(3))>frac
        fprintf('\n%2.0f%% Completed',frac*100);
        frac=frac+0.1;
    end
end

%Extraction
RefTAC  = zeros(nFrame,1);
RefMask  = zeros(sizePET(1:3));
v		= zeros(nFrame,1);
num_idx = 0; %COUNT THE NUMBER OF VOXELS IN THE REFERENCE

for z=1:sizePET(3)
    for x=1:sizePET(1)
        for y=1:sizePET(2)
            if (refMaskExtraction(x,y,z)>0) && (GRAYRATIO(x,y,z)>threshold)
                v(:) 		= dynPET(x,y,z,:); %using non-normalized PET for TAC extraction
                RefTAC		=  RefTAC+v;
                RefMask(x,y,z)	= 1; 
                num_idx=num_idx+1;
            end
        end
    end
end
fprintf('\nReference extraction complete\n');
RefTAC = RefTAC/num_idx;


end

