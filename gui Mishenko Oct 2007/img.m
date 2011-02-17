function img(X,Y,p)
% img(X,Y,p)
% a handy routine for displaying sum of images
% X and Y with weight coefficients (1-p) and p
% in grayscale

if(nargin<2)
    Y=0; p=0;
elseif(nargin<3)
    p=0.7;
end

h=figure;

if((nargin<2) | (~isequal(size(Y),size(X))))
    imagesc(X);
else
    imagesc(imlincomb(1-p,im2uint8(X),p,im2uint8(Y))); 
end

colormap gray;
