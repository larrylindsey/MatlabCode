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
if ~isfield(param, 'undImages')
    if isempty(param.trans)
        if isempty(param.imGrid)
            error('Both parameters trans and imGrid are empty');
        end
        
        im = getImage(param.imGrid);
        
        [rc_found, rc_grid, grid_model, scale] = ...
            extractDistortionTransform(im);
        
        exParam = getExtractedTransform();
        %exParam.norm = scale;
        exParam.order = param.transOrder;
        trans = getExtractedTransform(rc_found, rc_grid, grid_model,...
            exParam);
        
        transfile = sprintf('%s_trans.mat', param.prefix);
        save(transfile, 'trans', '-v7.3');
        param.trans = trans;
        clear trans;
    end
    
    % Transform images, then save to disk
    param.undImages = cell(size(param.images));
    
    for ii = 1:numel(param.images)
        im = getImage(param.images{ii});
        
        imtr = applyTransformImage(im, param.trans);
        
        param.undImages{ii} = sprintf('%s_und_%s',...
            param.prefix, param.images{ii});
        
        imwrite(imtr, param.undImages{ii});
    end
end

% Run Fiji
param.rawMatch = str2fun(runFiji(param.fijiExec, param.jarPath, ...
    param.images, [param.prefix '_raw']));
param.undMatch = str2fun(runFiji(param.fijiExec, param.jarPath, ...
    param.undImages, [param.prefix '_und']));

[rawA rawB] = param.rawMatch();
[undA undB] = param.undMatch();
[rawUndA rawUndB] = transformRaw(param.transform, rawA, rawB);

param.raw_rawError = computeError(rawA, rawB);
param.und_undError = computeError(undA, undB);
param.raw_undError = computeError(rawUndA, rawUndB);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

