function param = distortionApproxErrorAnalysis(param)

if nargin < 1
    param.trans = [];
    param.transOrder = 5;
    param.imGrid = '';
    param.images = {};
    param.prefix = '';
    param.fijiExec = '';
    param.jarPath = '~/code/java/';
    return
end

% Check inputs
errstr = '';

if isempty(param.prefix)
    errstr = ' :Prefix parameter cannot be empty';
end

if isempty(param.fijiExec)
    errstr = [errstr ' :Need Fiji executable'];
end

if ~isempty(errstr)
    error(errstr);
end

% If there is no image transform, make one
if ~isfield(param, 'undImages') || isempty(param.trans)
    if isempty(param.trans)
        if isempty(param.imGrid)
            error('Both parameters trans and imGrid are empty');
        end
        
        im = getImage(param.imGrid);
        
        [rc_found, rc_grid, grid_model, scale, db] = ...
            extractDistortionTransform(im);
        
        exParam = getExtractedTransform();
        %exParam.norm = scale;
        exParam.order = param.transOrder;
        exParam.useRansac = false;
        %exParam.affine = true;
        [trans junk junk terr] = getExtractedTransform(rc_found, rc_grid, grid_model,...
            exParam); %#ok
        
        transfile = sprintf('%s_trans.mat', param.prefix);        
        param.trans = trans;
        param.trans.scale = scale;
        param.trans.db = db;
        param.terr = terr;
        
        save(transfile, 'trans', '-v7.3');
        clear trans;
    end
    
    % Transform images, then save to disk
    param.undImages = transformImages(param.images, param.trans, ...
        param.prefix, param.trans.db);

end

% Run Fiji
param.rawMatch = str2func(runFiji(param.fijiExec, param.jarPath, ...
    param.images, [param.prefix '_raw']));
param.undMatch = str2func(runFiji(param.fijiExec, param.jarPath, ...
    param.undImages, [param.prefix '_und']));

[rawA rawB] = param.rawMatch();
[undA undB] = param.undMatch();
[rawUndA rawUndB] = transformRaw(param.trans, rawA, rawB);

param.rawPts.a = rawA;
param.rawPts.b = rawB;
param.undPts.a = undA;
param.undPts.b = undB;

param.err.raw_rawError = computeError(rawA, rawB);
param.err.und_undError = computeError(undA, undB);
param.err.raw_undError = computeError(rawUndA, rawUndB);

param.err.e_rr = cat(1, param.err.raw_rawError{:});
param.err.e_ru = cat(1, param.err.raw_undError{:});
param.err.e_uu = cat(1, param.err.und_undError{:});

[param.r.r_rr param.r.rxy_rr] = computeCorrcoef(rawA, rawB, param.images);
[param.r.r_uu param.r.rxy_uu] = ...
    computeCorrcoef(undA, undB, param.undImages);
[param.r.r_ru param.r.rxy_ru] = ...
    computeCorrcoef(rawUndA, rawUndB, param.undImages);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undImages = transformImages(images, trans, prefix, db)
undImages = cell(size(images));

msz = max(size(db.match)) * .5;

imlogi = true(size(getImage(images{1})));
imlogitr = not(applyTransformImage(imlogi, trans));
imlogitr([1 end],:) = true;
imlogitr(:,[1 end]) = true;
zeromask = bwdist(imlogitr) < msz;


parfor ii = 1:numel(images)
    im = getImage(images{ii});
    
    imtr = applyTransformImage(im, trans);
    
    %imtr(zeromask) = 0;
    
    [pstr fstr estr] = fileparts(images{ii});
    
    if isempty(pstr)
        pstr = '.';
    end
        
    undImages{ii} = sprintf('%s/%s_und_%s%s',...
        pstr, prefix, fstr, estr);
    
    imwrite(imtr, undImages{ii});
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r rxy] = computeCorrcoef(ptsA, ptsB, images)
r = cell(size(ptsA));
rxy = r;
for ii = 1:numel(ptsA)
    if ~isempty(ptsA{ii})
        [r{ii} rxy{ii}] = calculateAlignmentR(64, ptsA{ii}, ptsB{ii}, ...
           imread(images{ii}), imread(images{ii + 1}));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = computeError(ptsA, ptsB)
e = cell(size(ptsA));
for ii = 1:numel(e)
    if ~isempty(ptsA{ii})
        tr = regressionTransform(ptsA{ii}, ptsB{ii}, 1, @taylorMat);
        ptsAtr = doTransform(ptsA{ii}, tr);
        d = sqrt(sum((ptsAtr - ptsB{ii}).^2, 2));
        e{ii} = d;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A B] = transformRaw(tr, inA, inB)
A = cell(size(inA));
B = cell(size(inB));

scale = max(tr.scale);

for ii = 1:numel(A)
    if isempty(inA{ii})
        A{ii} = [];
        B{ii} = [];
    else
        scaleInA = applyScale(inA{ii}, scale);
        scaleInB = applyScale(inB{ii}, scale);
        
        A{ii} = reverseScale(doTransform(scaleInA, tr), scale);
        B{ii} = reverseScale(doTransform(scaleInB, tr), scale);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = applyScale(x, scale)
x = 2 * (x / scale) - 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = reverseScale(y, scale)
y = scale * (y + 1) / 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scriptfile = runFiji(fijiExec, jarpath, images, prefix)
% write the list file to disk
listfile = [prefix, '_list'];
scriptfile = [prefix, '_load'];
fid = fopen(listfile, 'w');
for ii = 1:numel(images)
    fprintf(fid, '%s\n', images{ii});
end
fclose(fid);

unixcmd = [fijiExec, ...
    ' --jarpath ', jarpath, ...
    ' --main-class SIFTMatchPointCollector ', ...
    listfile, ' ', scriptfile];
disp(['Executing command ' unixcmd])
tic;
unix(unixcmd);
toc;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = getImage(fileName)
% Read an image as a 2D single-precision floating point array
im = imread(fileName);
if size(im, 3) > 1
    im = rgb2gray(im);
end
im = im2single(im);
end

