function [eraw, eblock] = distortionTransformErrorBlock(trblock, cacheFile,...
    exclude)
% trblock - a block of transform structs, in the same shape as paramblock,
%           fit according to the correspoding parameters. eblock
%           corresponds, dimension-wise
% cachefile - sift cache file, as left by cacheFeatureMatch()
%
% eraw - the error output by analyzeSiftCache, using only an identity
%        transform
% eblock - the error output by analyzeSiftCache over transforms created
%           using paramblock.

if nargin < 5
    exclude = [];
end


if ischar(cacheFile)
    cacheFile = load(cacheFile);
end

eraw = analyzeSiftCacheAlignment(cacheFile, [], exclude);

eblock = analyzeSiftCacheAlignment(cacheFile, trblock(1), exclude);
eblock = repmat(eblock, size(trblock));


for ii = 2:numel(trblock)
    eblock(ii) = analyzeSiftCacheAlignment(cacheFile, trblock(ii), exclude);
end

