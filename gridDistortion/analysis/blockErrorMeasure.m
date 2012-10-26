function [emeasure, iorder] = blockErrorMeasure(eraw, eblock, fun, sortdir)
if nargin < 4
    sortdir = 'ascend';
end

if iscell(fun)
    efun = fun{1};
    xfun = fun{2};
else
    efun = fun;
    xfun = @(x) x.all;
end

emeasure = zeros(size(eblock));

for ii = 1:numel(eblock)
    emeasure(ii) = efun(xfun(eraw), xfun(eblock(ii)));
end

[~, iorder] = sort(emeasure(:), sortdir);
