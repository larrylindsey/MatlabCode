function [mask,p,s]=create_fold_mask3(I,th,area)

% Create a mask with background 0, and 1's where the image I has a region
% bigger than area and darker than th*max(I).
% If the fold doesn't cross completly the image, p is the coefficients of
% the line that cut the images in two subimages making the fold longer. 

% Example:
% [mask,p]=create_fold_mask2(Image,0.3,500);
% imshow(mask); 

thres=zeros(size(I));
thres(I/(max(I(:)))<th)=1;
mask=bwareaopen(thres,area);

%[Lregions,numregions]=bwlabel(imcomplement(mask),8);
[Lfolds,numfolds]=bwlabel(mask,8);
if max(mask(:))<1
    'There are no folds in this image.'
    p=0;
else
    p=zeros(numfolds,2);
    
    for i=1:numfolds
    [x,y]=find(Lfolds==i);
    [p(i,:),s(i)]=polyfit(x,y,1);
    end
    
end

    
        

% Marta Rivera-Alba
% 082908 Janelia Farm Research Campus. HHMI.