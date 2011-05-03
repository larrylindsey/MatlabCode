function [labelMask tdata] = selectLabelsByMask(label, mask, lt, ut)
% label - output from bwlabel or watershed
% mask - selection mask to match against label
% lt - lower threshold, used for nondirectional cover
% ut - upper threshold, used for directional cover

if numel(size(label)) ~= numel(size(mask)) ...
        || any(size(label) ~= size(mask))
    error('label and mask have different sizes, which is the poops');
end

if ~islogical(mask)
    error('mask must be logical');
end

L = bwlabel(mask);
n = max(L(:));

labelMask = false([size(mask) n]);
tdata = cell(1,n);

for il = 1:max(L(:))
    subMask = L == il;
    idL = double(unique(label(subMask)));
    idL(idL == 0) = [];
    innerScore = getInnerScore(label, subMask, idL);
    innerSelect = innerScore >= lt;
    
    outterScore = getOutterScore(label, subMask, idL, innerSelect);
    outterSelect = outterScore >= ut;
    
    tdata{il} = cat(2, idL(:), innerScore(:), innerSelect(:),...
        outterScore(:), outterSelect(:));
    
    labelSelect = cat(1, idL(innerSelect(:)), idL(outterSelect(:)));
    
    for ii = 1:numel(labelSelect)
        labelMask(:,:,il) = ...
            or(labelMask(:,:,il), label == labelSelect(ii));
    end
    
%     [labelMask tdata] = labelMaskHelper(label, L == il, lt, ut, labelMask,...
%         tdata);
end

labelMask = max(labelMask, [], 3);
tdata = cat(1, tdata{:});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = getInnerScore(label, mask, idL)
score = zeros(numel(idL), 1);
parfor ii = 1:numel(idL)
    currLabelMask = label == idL(ii);
    score(ii) = directionalCover(currLabelMask, mask);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outScore = getOutterScore(label, mask, idL, exclude)

outScore = zeros(numel(idL), 1);
idLi = idL(not(exclude));
idLe = idL(exclude);

score = zeros(numel(idLi), 1);

for ii = 1:numel(idLe)
    mask(label == idLe(ii)) = false;
end

parfor ii = 1:numel(idLi)
    currLabelMask = label == idLi(ii);
    score(ii) = nonDirectionalCover(currLabelMask, mask);
end

outScore(not(exclude)) = score;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [labelMask outdata] = labelMaskHelper(label, mask, lt, ut, ...
%     inLabelMask, indata)
%     
% lid = unique(label(mask));
% lid(lid == 0) = [];
% 
% % If a segment is approximately internal to the selection mask, add it to
% % the output.
% 
% innerSelect = false(size(lid));
% labelMask = false(size(mask));
% innerScore = zeros(size(lid));
% 
% for ii = 1:numel(lid)
%     currLabelMask = label == lid(ii);
%     innerScore(ii) = directionalCover(currLabelMask, mask);
%     if innerScore(ii) >= lt
%         innerSelect(ii) = true;
%         labelMask = or(labelMask, currLabelMask);
%     end
% end
% 
% innerLid = lid(innerSelect);
% innerScore = innerScore(innerSelect);
% innerLid = innerLid(:);
% innerScore = innerScore(:);
% 
% 
% outterSelect = false(size(lid));
% outterScore = zeros(size(lid));
% outterMask = mask;
% % Don't consider the contributions to the inner selections for this step.
% outterMask(labelMask) = false;
% 
% seq = 1:numel(lid);
% seq(innerSelect) = [];
% 
% for ii = seq
%     currLabelMask = label == lid(ii);
%     outterScore(ii) = nonDirectionalCover(currLabelMask, outterMask);
%     if outterScore(ii) >= ut
%         outterSelect(ii) = true;
%         labelMask = or(labelMask, currLabelMask);
%     end
% end
% 
% outterLid = lid(outterSelect);
% outterScore = outterScore(outterSelect);
% outterLid = outterLid(:);
% outterScore = outterScore(:);
% 
% labelMask = or(labelMask, inLabelMask);
% lid = cat(1, innerLid, outterLid);
% scores = cat(1, innerScore, outterScore);
% indices = cat(1, zeros(numel(innerScore), 1), ones(numel(outterScore), 1));
% 
% 
% outdata = cat(1, indata, cat(2, lid, scores, indices));
% end
