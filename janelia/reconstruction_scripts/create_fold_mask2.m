function [mask,p]=create_fold_mask2(I,th,area)

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

[L,num]=bwlabel(imcomplement(mask),8);
[L2,num2]=bwlabel(mask,8);
 
if max(mask(:))<1
    'There are no folds in this image.'
    p=0;
elseif (num==1 && num2==1)
    
    [x,y]=find(mask>0);
    p=polyfit(x,y,1);
end

if num2>1
    'There are more than one fold in this image'
    
end


    
        

% Marta Rivera-Alba
% 082908 Janelia Farm Research Campus. HHMI.