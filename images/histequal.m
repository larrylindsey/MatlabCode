function imout = my_histequal(im)
if ~isinteger(im)
    error('This function only works on integer images at present');
end

if size(im, 3) > 1
    error('This function only works on greyscale images at present');
end

[hh] = myhist(im);

hhnorm = hh / numel(im);
HH = cumsum(hhnorm);

%Here, we trade memory-efficiency for speed.
imcopy = repmat(im, [1 1 256]);
rpsz = [size(im) 1];

checkmat = uint8(zeros(1, 1, 256));
checkmat(:) = 0:255;
checkmat = repmat(checkmat, rpsz);

valmat = zeros(1, 1, 256);
valmat(:) = HH;
valmat = repmat(valmat, rpsz);

checkmat = double(imcopy == checkmat);
valmat = sum(valmat .* checkmat, 3);

imout = im2uint8(valmat);

end