function [stackOut matchData score] = autoSegmentReconstruct(secdoc, ...
    cname, index, cstr)

if nargin < 4
    cstr = defaultControl;
    if nargin == 0
        stackOut = cstr;
        return;
    else
        [x y] = reconstructDomainBounds(secdoc);
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
nextMask = maskOrig;
nextData = dataOrig;
nextScore = -1;

doRepartition = false;


while any(nextMask(:)) || doRepartition
    
    if ~doRepartition
        
        stackOut = cat(3, stackOut, nextMask);
        matchData = cat(2, matchData, nextData);
        score = cat(2, score, nextScore);
        currMask = nextMask;
        
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
    
    if nextScore < cstr.scoreT        
        return;
        if doRepartition
            fprintf('Repartitioning %d led to no further segmentation\n',...
                index);
            return;
        end
        disp('Score did not meet requirement. Resegmenting');
        doRepartition = true;
        
        label = nCutPartition(label, currMask, nextData, index);
    else
        doRepartition = false;
    end
end

stackOut = cat(3, stackOut, nextMask);
matchData = cat(2, matchData, nextData);
score = cat(2, score, nextScore);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partLabel = nCutPartition(label, mask, data, index)
partLabel = label;
% This usually should happen because the mask is much smaller than an
% under-segmented component in the label that it would otherwise correspond
% to.
% Here, find that segment by figuring out which component the mask fits
% best within.

n = 4;

idL = data(:,1);
matchScore = zeros(size(idL));

parfor il = 1:numel(idL)
    matchScore(il) = directionalCover(mask, label == idL(il));
end

[junk iMatch] = max(matchScore);

segMask = label == idL(iMatch);

classIm = readClassIm(index);
classIm = sizeImage(classIm, size(segMask));

[W junk cropMask offset] = nCutMaskedImage(segMask, classIm, 4, [1 1]);
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

repMask = false(size(label));
repMask(offset(1) + (1:size(cropMask,1)), offset(2) + (1:size(cropMask,2))) ...
    = cropMask;

partLabel(repMask) = newLabel(cropMask);
partLabel(repMask) = partLabel(repMask) + max(partLabel(:));

partLabel(partLabel == max(partLabel(:))) = idL(iMatch);


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