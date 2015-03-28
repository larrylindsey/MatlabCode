function im = replaceColor(im, c1, c2)

if ~all(size(c1) == size(c2))
    error('Colors must have the same shape');
end

im = im2double(im);

k = size(c1, 1);

im_r = im(:,:,1);
im_g = im(:,:,2);
im_b = im(:,:,3);

for i_k = 1:k
    c_in = c1(i_k,:);
    c_out = c2(i_k,:);
    mask = getMask(im, c_in); 
    im_r(mask) = c_out(1);
    im_g(mask) = c_out(2);
    im_b(mask) = c_out(3);
end

im = cat(3, im_r, im_g, im_b);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = getMask(im, c)
% Creates a mask over im for a given color

% Start out everywhere true.
mask = true(size(im,1), size(im,2));
for ic = 1:numel(c)
    % Set the mask to false if the image's channel
    % doesn't match the corresponding component of
    % the color. This works for RGB or grayscale.
    mask = and(mask, abs(im(:,:,ic) - c(ic)) < (16 / 256));
end


end