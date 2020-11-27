function IDIF = computeIDIF(PET,delta,timeToPeak, BRAIN);
nVoxel = 40; %BASED ON TOMASI ORIGINAL CODE
%DEFINE THE TIME INTERVAL and SUM IMAGE
endFrame=cumtrapz(delta);
ind=find(endFrame<timeToPeak); 
sumPET = sum(PET(:,:,:,ind),4);
sumPET = sumPET.*BRAIN; 
%FIND THE N VOXELS WITH THE HIGHEST ACTIVITY
indIDIF=zeros(nVoxel,1);
for n = 1:nVoxel
    [~,indIDIF(n)] = max(sumPET(:));
    sumPET(indIDIF(n))=-1;
end
%DEFINE THE IDIF
IDIF=zeros(length(delta),1);
for f=1:length(IDIF)
    temp=squeeze(PET(:,:,:,f));
    IDIF(f)=mean(temp(indIDIF));
end

