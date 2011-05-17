function imc = normcorr(im, kern)

if ~isfloat(im)
    im = im2single(im);
end

if ~isfloat(kern)
    kern = im2single(kern);
end

avg = ones(size(kern));

%kern = kern - mean(kern(:));

kern = kern / norm(kern(:));

%imavg = fftfilter2(im, avg);
%im = im - imavg;

imnorm = sqrt(fftfilter2(im.*im, avg));

imkf = fftfilter2(im, kern);

imc = imkf ./ imnorm;