function gridTem = getRealTEMGrid(imTem, imStem, gridTem, ramStem, gridStem)
% function outGridTem = 
%                      getRealTEMGrid(imTem, imStem, gridTem, ramStem, gridStem)
%   imTem - a TEM calibration image
%   imStem - an STEM calibration image, in which imTem is inset.
%   gridTem - a calibration struct for imTem, containing the following
%             fields:
%       rc_found - the row-column location of alignment features, typically
%                  the corners between cross-grating cells.
%       rc_grid - the square-grid locations corresponding to rc_found
%       grid_model - an approximated 2D vector basis over rc_found
%           These fields should be as returned by
%           extractDistortionTransform
%   ramStem - a calibration struct over a known-rigid image taken in the
%             Stem. As of this writing, this will likely be over a
%             backscatter image of a RAM chip.
%   gridStem - a calibration struct over imStem.
%
%   outGridTem - a calibration struct like gridTem, with the following
%                differences:
%       - rc_found will be decimated so that every location, as found in
%         imTem, corresponds to the identical location found in imStem.
%       - a field, rc_correct, is added. rc_correct represents the same
%         locations as in rc_found as they are detected in imStem, after
%         correction for distortion. In other words, these locations are
%         where we *should* see the discovered features, if there were no
%         distortion in the TEM.

imStem = verifyImage(imStem);
imTem = verifyImage(imTem);
gridTem = verifyStruct(gridTem);
gridStem = verifyStruct(gridStem);
ramStem = verifyStruct(ramStem);

[imStem, r] = doScale(imStem);
if r ~= 1
    imTem = imresize(imTem, r);
end

if ischar(gridTem)
    gridTem = load(gridTem, 'rc_found', 'rc_grid', 'grid_model');
end

if ischar(gridStem)
    gridStem = load(gridStem, 'rc_found', 'rc_grid', 'grid_model');
end



[rc_tem_match, rc_stem_match] = getMatchSet(imTem, imStem);
trRam = getRamTransform(ramStem);
[gridTem.rc_found, gridTem.rc_grid, gridTem.rc_correct] = ...
    getRCReal(trRam, rc_tem_match, rc_stem_match, gridTem, gridStem);

gridTem = affineAlignRCCorrect(gridTem);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = verifyImage(im)
if ischar(im)
    im = imread(im);
end

im = im2single(im);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = verifyStruct(str)
if ~isstruct(str)
    str = load(str, 'rc_found', 'rc_grid', 'grid_model');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [im, r] = doScale(im)
MAXSIZE = 4096;
f = round(max(size(im)) / MAXSIZE);
r = 1 / f;
im = imresize(im, r);
fprintf('Image scaled to size %d by %d\n', size(im,1), size(im,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rc_tem, rc_stem] = getMatchSet(imTem, imStem)
% extract SIFT features from both images.
fprintf('Extracting SIFT features from STEM...\n');
[f_stem, d_stem] = thirdpartyfunction('vl_sift', imStem);
fprintf('Extracting SIFT features from TEM...\n');
[f_tem, d_tem] = thirdpartyfunction('vl_sift', imTem);

fprintf('Matching features...\n');
% find putative matches
m_s_t = thirdpartyfunction('vl_ubcmatch', d_stem, d_tem);
% since we're matching an inset image, we can perform a first-pass
% rejection based on the match density.
[f_stem, f_tem] = sparseMatchReject(f_stem, f_tem, 2, m_s_t);

% Scale the matches to [-1 1] x [-1 1]. This is done to be size invariant,
% and also because extractDistortionTransform returns rc_found in the same
% space.
rc_stem = f_stem([2 1],:)' * 2 / max(size(imStem)) - 1;
rc_tem = f_tem([2 1],:)' * 2 / max(size(imTem)) - 1;

rparam.n = 7;
rparam.metric = [];
rparam.maxError = .1;
rparam.minInliers = 32;
rparam.maxIter = 1024;

% trInlay will be used to transform the points rc_found, gotten by running
% extractDistortionTransform over the TEM image, into the space of the STEM
% image.
[~, sel] = ransacRegressionTransform(rparam, rc_tem, rc_stem, 1);
rc_tem = rc_tem(sel,:);
rc_stem = rc_stem(sel,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trRam = getRamTransform(ramStem)
control = getExtractedTransform;
control.order = 5;
control.xAngle = pi / 3;
control.useRansac = false;
trRamStr = getExtractedTransform(ramStem.rc_found, ramStem.rc_grid, ...
    ramStem.grid_model, control);
trRam = trRamStr.affine.tr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rcTem, rcTemGrid, rcStem] = getRCReal(trRam, rc_tem_match, rc_stem_match, ...
    gridTem, gridStem)
rc_stem_match = doTransform(rc_stem_match, trRam);
tr_tem_stem = regressionTransform(rc_tem_match, rc_stem_match, 1);
tr_stem_tem = regressionTransform(rc_stem_match, rc_tem_match, 1);

[rcTem, rcTemGrid, rcStem] = selectMatchingGridPoints(tr_tem_stem, gridTem, gridStem);

rcStem = doTransform(rcStem, tr_stem_tem);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rc_found_tem rc_grid_tem rc_found_stem] = ...
    selectMatchingGridPoints(tr_tem_stem, gridTem, gridStem)
bbox_buffer = max(abs(gridTem.grid_model(:)));
tem_bbox = [-1 -1; -1 1; 1 1; 1 -1];
tem_bbox(tem_bbox > 0) = tem_bbox(tem_bbox > 0) + bbox_buffer;
tem_bbox(tem_bbox < 0) = tem_bbox(tem_bbox < 0) - bbox_buffer;

stem_bbox = doTransform(tem_bbox, tr_tem_stem);
% Use a bounding box to reduce the search size for matching the individual...
% grid points TEM vs STEM.
stem_coarse_sel = inpolygon(gridStem.rc_found(:,1), gridStem.rc_found(:,2), ...
    stem_bbox(:,1), stem_bbox(:,2));
rc_found_stem = gridStem.rc_found(stem_coarse_sel,:);
rc_found_tem = doTransform(gridTem.rc_found, tr_tem_stem);

% Do the actual selection by nearest neighbor
dd_thresh = max(sqrt(sum(gridStem.grid_model.^2, 2)));
dd = dist2(rc_found_tem, rc_found_stem);
[ddsort, iddsort] = sort(dd, 2, 'ascend');
ddsort = ddsort(:,1);

ok_sel = sqrt(ddsort) * 2 < dd_thresh;

% assume our grids are good enough to do a one to one mapping
stem_sel = iddsort(ok_sel ,1);
if numel(stem_sel) ~= numel(unique(stem_sel))
    warning('Oh, shit'); %#ok<WNTAG>
    keyboard;
end

rc_found_tem = gridTem.rc_found(ok_sel,:);
rc_grid_tem = gridTem.rc_grid(ok_sel,:);
rc_found_stem = rc_found_stem(stem_sel,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gridTem = affineAlignRCCorrect(gridTem)
rc_correct = gridTem.rc_correct;
rc_found = gridTem.rc_found;
trAff = regressionTransform(rc_correct, rc_found, 1);
rc_correct_aff = doTransform(rc_correct, trAff);
gridTem.rc_correct_aff = rc_correct_aff;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
