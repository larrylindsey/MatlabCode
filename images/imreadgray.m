function im = imreadgray(file)

im = imread(file);
if size(im, 3) > 1
    im = rgb2gray(im);
end
