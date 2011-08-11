function [im_out sz_out seams_out] = reduceHeight(im, varargin)
%function [im_out sz_out seams_out] = reduceHeight(im, [numPixels, [efunc]])
%Reduces the height of an image by a certain number of pixels, using the
%seam carving algorithm.
%
% im - the input image.  May be integer or double.
% numPixels - the number of rows to remove from the image.  Default: 1
% efunc - a function handle to a function that takes a double image with
%         intensity range [0 1], and outputs an energy map for that image.
%         The default is to use sumOfDerivativesEnergy
%
% im_out - the output image, with the specified number of rows removed.
%          This image will be double-valued with intensity range [0 1]
% sz_out - [2 numPixels] is a pair of vectors representing the size history
%          of the image.  sz_out(1,:) is columns, while sz_out(2,:) is
%          rows.
% seams_out - a cell array containing the seams removed.  The seams are
%             represented as vectors, where the i'th element indexes the
%             row that was removed at column i.  The seams are stored in
%             the cell array in the order in which they were removed.  Note
%             that for each successive seam, the image is 1 pixel smaller
%             in height.
%
% This function is a wrapper around reduceWidth.  It feeds a
% transpose-image to that function, then adapts the output data to a height
% reduction operation.
%
%See also sumOfDerivativesEnergy reduceWidth

if nargout > 2
    [im_out sz_out seams_out] = ...
        reduceWidth(permute(im, [2 1 3]), varargin{:});
else
    [im_out sz_out] = reduceWidth(permute(im, [2 1 3]), varargin{:});
end

im_out = permute(im_out, [2 1 3]);
sz_out = sz_out(:, [2 1]);