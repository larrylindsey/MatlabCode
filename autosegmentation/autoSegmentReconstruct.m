function [stackOut matchData score indexout] = autoSegmentReconstruct(secdoc, ...
    cname, index, cstr)
[x y] = reconstructDomainBounds(secdoc);
if nargin < 4
    cstr = defaultControl;
    if nargin == 0
        stackOut = cstr;
        return;
    else
        
        cstr.x = x;
        cstr.y = y;
        cstr.eT = .5;
        cstr.iT = .8;
        cstr.scoreT = .75;
    end
end


% Get first match mask
annoMask = getRawImageContourSlice(secdoc, index, cname);
trOrig = secdocTrans(secdoc, index);
labelOrig = readLabelIm(index);
cstr.imsz = size(labelOrig);

annoMaskTR = applyTransformImage(annoMask, trOrig, x, y);
labelOrig = sizeImage(labelOrig, size(annoMaskTR));


[maskOrig dataOrig] = ...
    selectLabelsByMask(labelOrig, annoMaskTR, cstr.iT, cstr.eT);


stackOut = false(0);
matchData = {};
score = [];
indexout = [];
nextMask = maskOrig;
nextData = dataOrig;
nextScore = -1;

doRepartition = false;

finish = false;

while (any(nextMask(:)) || doRepartition) && ~finish
    
    if ~doRepartition
        
        stackOut = cat(3, stackOut, nextMask);
        matchData = cat(2, matchData, nextData);
        score = cat(2, score, nextScore);
        currMask = nextMask;
        indexout = cat(2, indexout, index);
        
        index = index - 1;
        
        fprintf('Now segmenting index %d\n', index);
        
        label = readLabelIm(index);
        label = sizeImage(label, size(currMask));    

    else
        fprintf('Re segmenting index %d\n', index);
    end

    
    [nextMask nextData] = selectLabelsByMask(label, currMask, cstr.iT,...
        cstr.eT);
    
    nextScore = nonDirectionalCover(nextMask, currMask);
    
    fprintf('Match score: %g\n', nextScore);
    
    if nextScore < cstr.scoreT  
%         keyboard
        if doRepartition
            fprintf('Repartitioning %d led to no further segmentation\n',...
                index);
            finish = true;
        else
            disp('Score did not meet requirement. Re partitioning');
            doRepartition = true;
            
            label = nCutPartition(label, currMask, nextData, index);
            moveLabelIm(index);
            writeLabelIm(label, index);
        end
    else
        doRepartition = false;
    end
end

stackOut = cat(3, stackOut, nextMask);
matchData = cat(2, matchData, nextData);
score = cat(2, score, nextScore);
indexout = cat(2, indexout, index);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partLabel = nCutPartition(label, mask, data, index)
partLabel = label;
% This usually should happen because the mask is much smaller than an
% under-segmented component in the label that it would otherwise correspond
% to.
% Here, find that segment by figuring out which component the mask fits
% best within.

idL = data(:,1);
matchScore = zeros(size(idL));

parfor il = 1:numel(idL)
    matchScore(il) = directionalCover(mask, label == idL(il));
end

[junk iMatch] = max(matchScore);

segMask = label == idL(iMatch);

classIm = readClassIm(index);
classIm = sizeImage(classIm, size(segMask));


cutMask = and(segMask, classIm > .75);

cutMaskSplit = splitCutMask(cutMask);

partLabel(segMask) = 0;

for ic = 1:size(cutMaskSplit,3)
    partLabel = nCutHelper(partLabel, cutMaskSplit(:,:,ic), classIm);
end

partLabel = fillNN(partLabel, segMask, cutMask);

partLabel(partLabel == max(partLabel(:))) = idL(iMatch);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maskSplit = splitCutMask(mask)
[maskLabel nl] = bwlabel(mask,4);

if nl > 2
    lhist = zeros(1, nl);
    nzMaskLabel = maskLabel(maskLabel > 0);
    for ii = 1:numel(nzMaskLabel)
        k = nzMaskLabel(ii);
        lhist(k) = lhist(k) + 1;
    end
    
    lhistSort = sort(lhist, 'descend');
    
    cullSel = find(lhist < mean(lhistSort(1:2)) / 16);
    
    for ii = 1:numel(cullSel)
        maskLabel(maskLabel == cullSel(ii)) = 0;
    end
    
    [maskLabel nl] = bwlabel(maskLabel > 0, 4);
elseif nl == 2
    lhist = zeros(1, 2);
    nzMaskLabel = maskLabel(maskLabel > 0);
    for ii = 1:numel(nzMaskLabel)
        k = nzMaskLabel(ii);
        lhist(k) = lhist(k) + 1;
    end
    
    if max(lhist) > min(lhist) * 16
        ii = find(lhist == min(lhist));
    end
    
   maskLabel(maskLabel == ii) = 0;
   
   maskLabel(maskLabel > 0) = 1;
   
   nl = 1;
    
end

maskSplit = false([size(mask) nl]);

for ii = 1:nl
    maskSplit(:,:,ii) = maskLabel == ii;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partLabel = nCutHelper(partLabel, mask, classIm)
n = 4;
[W junk cropMask offset] = nCutMaskedImage(mask, classIm, 4, [1 1]);
offset = offset - 1;

[e v d] = doNCuts(W, n);

segValue = zeros([size(cropMask) n]);

for ii = 1:n
    eSeg = zeros(size(cropMask));
    eSeg(cropMask) = d * v(:,ii);
    segValue(:,:,ii) = eSeg;
end

%segValue = [];

segSign = prod(sign(segValue(:,:,e ~= min(e))), 3);

posSeg = segSign > 0;
negSeg = segSign < 0;

posSeg = imerode(posSeg, strel('disk', 3));
negSeg = imerode(negSeg, strel('disk', 3));

posLabel = bwlabel(posSeg, 4);
negLabel = bwlabel(negSeg, 4);

negLabel(negLabel > 0) = negLabel(negLabel > 0) + max(posLabel(:));

newLabel = negLabel + posLabel;
newLabel = imdilate(newLabel, strel('disk', 1));

repMask = false(size(partLabel));
repMask(offset(1) + (1:size(cropMask,1)), offset(2) + (1:size(cropMask,2))) ...
    = cropMask;


newLabel(newLabel > 0) = newLabel(newLabel > 0) + ...
    double(max(partLabel(:)));
partLabel(repMask) = newLabel(cropMask);
%partLabel(repMask) = partLabel(repMask) + max(partLabel(:));


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function classIm = readClassIm(index)
classIm = im2double(readCacheIm('class', index));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labelIm = readLabelIm(index)
labelIm = readCacheIm('label', index);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cacheIm = readCacheIm(prefix, index)
cacheIm = imread(sprintf('seg_cache/%s_%04d.png', prefix, index));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function moveLabelIm(index)
infile = sprintf('seg_cache/label_%04d.png', index);
outIndex = 0;
outfile = sprintf('seg_cache/label_%04d.png.%d', outIndex);
while exist(outfile, 'file') > 0
    outIndex = outIndex + 1;
    outfile = sprintf('seg_cache/label_%04d.png.%d', outIndex);
end
movefile(infile, outfile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeLabelIm(label, index)
imwrite(label, sprintf('seg_cache/label_%04d.png', index));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cstr = defaultControl
cstr.x = [];
cstr.y = [];
cstr.eT = .5;
cstr.iT = .8;
cstr.scoreT = .75;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = secdocTrans(secdoc, index)
if nargin > 1
    section = secdoc([secdoc.index] == index).section;
else
    section = secdoc;
end
imIndex = section.transImageIndex;
t = section.Transform(imIndex);

end