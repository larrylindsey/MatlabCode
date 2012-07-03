function [eraw, eblock, trblock] = distortionTransformErrorBlock(...
    rc_found, rc_match, paramblock, cacheFile)
% rc_found - rc_found as returned by extractGridDistortion
% rc_match - the locations to match rc_found. We want tr: rc_found -> rc_match
% paramblock - an array of struct with the following fields:
%              order - the order parameter passed into fitTransform
%              fun - a function handle, as passed into fitTransform
% cachefile - sift cache file, as left by cacheFeatureMatch()
%
% eraw - the error output by analyzeSiftCache, using only an identity
%        transform
% eblock - the error output by analyzeSiftCache over transforms created
%           using paramblock.
% trblock - a block of transform structs, in the same shape as paramblock,
%           fit according to the correspoding parameters. eblock
%           corresponds, dimension-wise


eraw = analyzeSiftCacheAlignment(cacheFile);

trblock = fitTransform(rc_found, rc_match, paramblock(1).order,...
    paramblock(1).fun);
eblock = analyzeSiftCacheAlignment(cacheFile, trblock);

trblock = repmat(trblock, size(paramblock));
eblock = repmat(eblock, size(paramblock));


parfor ii = 2:numel(paramblock)
    trblock(ii) = fitTransform(rc_found, rc_match, paramblock(ii).order,...
        paramblock(ii).fun);
    eblock(ii) = analyzeSiftCacheAlignment(cacheFile, trblock(ii));
end

