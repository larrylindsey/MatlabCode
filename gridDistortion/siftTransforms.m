function trArray = siftTransforms(files, tr)
% function trArray = siftTransforms(files, [tr])
% Collects sift-calculated transforms for matches between the files in the
% cell-array passed int
%
% If specified, the transform tr is applied to each image before sift
% analysis.

if isempty(which('vl_sift'))
    error('Please add vl_sift fo the path');
end

if nargin < 2
    tr = [];
end

rparam = ransacRegressionTransform;
rparam.n = 10;
rparam.maxError = 100;
rparam.minInliers = 12;
rparam.maxIter = 1000;

trArray = repmat(struct, [numel(files) - 1, 1]);

f = cell(numel(files));
d = f;

parfor ii = 1:numel(files)
    fprintf('Reading %s\n', files{ii});
    im = getImage(files{ii}, tr);    
    fprintf('Done reading, doing the SIFT on %s\n', files{ii});
    [f{ii} d{ii}] = vl_sift(im);
    fprintf('Done sifting %s\n', files{ii});
end

save -v7.3 siftcache f d

fcurr = f(2:end);
flast = f(1:(end - 1));
dcurr = d(2:end);
dlast = d(1:(end - 1));

%imcurr = getImage(files{2});
for ii = 1:(numel(files) - 1)
    fprintf('Finding Matches\n');
    matches = vl_ubcmatch(dcurr{ii}, dlast{ii});
    fprintf('Done.\n');
    
    matchLocLast = flast{ii}(1:2, matches(2,:))';
    matchLocCurr = fcurr{ii}(1:2, matches(1,:))';
    
    matchDist = sqrt(sum((matchLocLast - matchLocCurr).^2, 2));
    
    matchDistSort = sort(matchDist);
    th = matchDistSort(round(end / 8)); %first octile distance
    
    sel = matchDist < th;
    
    ptsLast = matchLocLast(sel,:);
    ptsCurr = matchLocCurr(sel,:);
    
    [trSift trSiftSel] =...
        ransacRegressionTransform(rparam, ptsLast, ptsCurr, 1);
    
    trSift.fromPts = ptsLast;
    trSift.toPts = ptsCurr;
    trSift.iFrom = ii - 1;
    trSift.iTo = ii;
    trSift.trSiftSel = trSiftSel;
    
    flast = fcurr;
    dlast = dcurr;
    
    trArray(ii) = trSift;
end



end

function im = getImage(path, tr)
im = imread(path);
if size(im, 3) > 1
    im = rgb2gray(im);
end

im = im2single(im);

if ~isempty(tr)
    im = applyTransformImage(im, tr);
end

end