function [image_out] = remove_seam(image_in, seamc)
%[image_out] = remove_seam(image_in, seamc)
%
%Removes a column-seam from the input image.
%
% image_in - the image from which a seam is to be removed.
% seamc - a vector that has as many elements as the image has rows.  The
%         element at the i'th position should correspond to the i'th row of
%         the image, and should contain a value that indexes the columnar
%         position of the pixel to be removed at that row.

imsz = size(image_in);
newimsz = imsz;

%The new image will have one fewer columns
newimsz(2) = newimsz(2) - 1;

%Basic find-type selector
basic_csel = 1:imsz(2);

%Pre-allocate the output image.
image_out = zeros(newimsz);

for r = 1:imsz(1)
    %Remove the index of the column in the seam at this row.
    csel = basic_csel(logical(basic_csel ~= seamc(r)));
    
    %Push the new row (with seam intersection removed) into the new image.
    image_out(r, :, :) = image_in(r,csel,:);
end

end