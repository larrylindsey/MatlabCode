function [imout id] = read_sop_image(path, sz)

idash = find(path == '_');
idot = find(path == '.');


id = str2num(path((idash(1) + 1):(idash(2) - 1)));

x = str2num(path((idash(2) + 1):(idash(3) - 1))) + 1;
y = str2num(path((idash(3) + 1):(idot - 1))) + 1;

imout = zeros(sz);

im = im2double(imread(path));

imout(y:(y + size(im,1) - 1), x:(x + size(im,2)) - 1) = im;

end
