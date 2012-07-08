function cacheSIFTFeatures(files, cachefile)

if isstruct(files)
    files = {files.name};
end

if nargin < 2
    iSlash = find(files{1} == '/', 1, 'last');
    if isempty(iSlash)
        iSlash = 0;
    end
    
    cachefile = [files{1}(1:iSlash) 'siftCache.mat'];
end

fprintf('Saving sift features and locations to %s\n', cachefile);

f = cell(1, numel(files));
d = f;
sz = zeros(numel(files), 2);

parfor ii = 1:numel(files)
    imfile = files{ii};   
    im = im2single(imread(imfile));
    if size(im,3) > 1
        im = im(:,:,1);
    end
    
    sz(ii,:) = size(im); %#ok<PFOUS>
    
    [f{ii} d{ii}] = thirdpartyfunction('vl_sift', im, ...
        'FirstOctave', 2); %#ok<PFOUS> 
end

% fprintf('f has size %d, and d has size %d\n', numel(f), numel(d));

save(cachefile, '-v7.3', 'f', 'd', 'sz', 'files');

