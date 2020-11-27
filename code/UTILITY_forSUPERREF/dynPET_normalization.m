function [ normDYN ] = dynPET_normalization(DYN,map)
%the function normalizes the PET dynamic data respect to the voxel define in map

DIM=size(DYN);
normDYN=zeros(size(DYN));
ind=find(map>0);

for f=1:1:DIM(4);
 temp=DYN(:,:,:,f);
 mF=mean(temp(ind));
 stdF=std(temp(ind));
 nF=(temp-mF)./stdF;
 normDYN(:,:,:,f)=nF;
 clear temp nF
end



