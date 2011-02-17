function [r t] = getFijiHoughPeaks(imvote, im, n)

pk = findPeaks2(imvote);

pkval = imvote(pk);
pkval = unique(pkval);
pkval = sort(pkval, 'descend');


%ord = ord(sel);

r = [];
t = [];

if n > numel(pkval)
    n = numel(pkval);
end

for i_n = 1:n
    [tt, rr] = find(imvote == pkval(i_n));
    r = cat(1, r, rr);
    t = cat(1, t, tt);
    
end

t = t - 91;

r = (r - 1 - (size(imvote, 2) / 2)) * sqrt(sum(size(im).^2)) * 2 / size(imvote, 2);