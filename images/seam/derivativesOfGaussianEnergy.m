function efunc = derivativesOfGaussianEnergy(image_in_dbl, h, s)
%derivativesOfGaussianEnergy(image_in_dbl, [h, [s]])
%
%Calculates an energy function over the input image.  Similar to
%sumOfDerivativesEnergy, but with an additional Gaussian filter.
%
% image_in_dbl - the input image.  Expected to be double-valued.
% h - the size of the gaussian filter to use.  Default: 12
% s - the standard deviation of the gaussian filter. Default: .5
%
% efunc - the energy computation output.
%
% See also sumOfDerivativesEnergy

if nargin < 3
    if nargin < 2
        h = 12;
    end
    s = .5;
end


image_in_gray = mean(image_in_dbl,3);
d_dx = [-1 1];
d_dy = [-1 ; 1];
gaussian = fspecial('gaussian', h, s);
gauss_dx = imfilter(gaussian, d_dx, 'same', 'replicate');
gauss_dy = imfilter(gaussian, d_dy, 'same', 'replicate');

efunc = zeros(size(image_in_gray));
efunc = efunc + abs(imfilter(image_in_gray, gauss_dx, 'same', 'symmetric'));
efunc = efunc + abs(imfilter(image_in_gray, gauss_dy, 'same', 'symmetric'));
