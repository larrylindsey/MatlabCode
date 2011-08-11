function efunc = sumOfDerivativesEnergy(image_in_dbl)
%efunc = sumOfDerivativesEnergy(image_in_dbl)
%
%Calculates an energy function over the input image.  The energy function
%is determined by the sum of the absolute vaues of the row-direction and
%column-direction intensity derivates, as determined by filtering with
%[-1 1] and [-1 ; 1].  These derivatives are calculated from a grayscale
%representation.
%
% image_in_dbl - the input image.  This is expected to be double-valued.
% 
% efunc - the energy computation output.  This is a two-dimensional matrix
%         with the same number of rows and colums as the input image.
%
% See also derivativesOfGaussianEnergy

image_in_gray = mean(image_in_dbl,3);
d_dx = [-1 1];
d_dy = [-1 ; 1];
efunc = zeros(size(image_in_gray));
efunc = efunc + abs(imfilter(image_in_gray, d_dx, 'same', 'symmetric'));
efunc = efunc + abs(imfilter(image_in_gray, d_dy, 'same', 'symmetric'));
