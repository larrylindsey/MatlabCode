function [samples labeledPeaks] = histogramToFenceStuff(histogram)

sz = size(histogram);

if sz(1) > sz(2)
    histogram = histogram';
    sz = size(histogram);
end

gg = fspecial('gaussian', 11, 5);

%L = bwlabel(histogram > 0);

% Remove jitter.
histogram = imfilter(histogram, gg, 'same');
L = watershed(histogram);
pk = findPeaks(histogram);

samples = zeros(numel(pk), 2);
samples(:, 1) = pk;
samples(:, 2) = L(pk);

labeledPeaks = zeros(numel(histogram), 2);
labeledPeaks(:, 1) = L;
labeledPeaks(:, 2) = histogram;
