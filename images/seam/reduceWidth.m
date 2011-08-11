function [im_out sz_out seams_out] = reduceWidth(im, numPixels, efunc, ...
    costmap)
%function [im_out sz_out seams_out] = 
%                          reduceWidth(im, [numPixels, [efunc, [costmap]]])
%Reduces the width of an image by a certain number of pixels, using the
%seam carving algorithm.
%
% im - the input image.  May be integer or double.
% numPixels - the number of columns to remove from the image.  Default: 1
% efunc - a function handle to a function that takes a double image with
%         intensity range [0 1], and outputs an energy map for that image.
%         The default is to use sumOfDerivativesEnergy
% costmap - if a costmap is passed in, it will be used for the first seam
%           removed, rather than calculating a new one.
%
% im_out - the output image, with the specified number of columns removed.
%          This image will be double-valued with intensity range [0 1]
% sz_out - [2 numPixels] is a pair of vectors representing the size history
%          of the image.  sz_out(1,:) is columns, while sz_out(2,:) is
%          rows.
% seams_out - a cell array containing the seams removed.  The seams are
%             represented as vectors, where the i'th element indexes the
%             column that was removed at row i.  The seams are stored in
%             the cell array in the order in which they were removed.  Note
%             that for each successive seam, the image is 1 pixel smaller
%             in width.
%See also sumOfDerivativesEnergy reduceHeight



%Set defaults
if nargin < 3
    efunc = @defaultEnergy;
    if nargin < 2
        numPixels = 1;
    end
end

%If the image is integer-valued, cast to double and scale down.
if isinteger(im)
    im_out= double(im) / 255;
else
    im_out = im;
end

%Pre-allocate memory
sz_out = zeros(numPixels, 2);
%Only allocate the seams if they were requested.  This could potentially be
%rather memory-intensive.
if nargout > 2
    seams_out = cell(1, numPixels);
end

%Main loop.  Hopefully self-explanatory.
for ii = 1:numPixels
    energy = efunc(im_out);
    %Special case in which cost map has been calculated beforehand
    if ii == 1 && nargin > 3
        seam = calculate_seam(energy, costmap);
    else
        seam = calculate_seam(energy);
    end

    if nargout > 2
        seams_out{ii} = seam;
    end
    im_out = remove_seam(im_out, seam);
    
    %Yeah, I could just write this as a vector, but this might also be
    %useful for debug purposes.
    sz = size(im_out);
    sz_out(ii, :) = sz(1:2);
end

end