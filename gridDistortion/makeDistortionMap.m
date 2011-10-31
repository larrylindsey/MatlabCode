function [dmag dcolor RC RCtr] = makeDistortionMap(trans, n, lim)

if nargin < 3
    lim = [-1 1];
end


[R C] = meshgrid(linspace(lim(1), lim(2), n), linspace(lim(1), lim(2), n));
RC = cat(2, C(:), R(:));
RCtr = doTransform(RC, trans);


dmag = reshape(sqrt(sum((RC - RCtr).^2, 2)), size(R));

dcolor = vectorMapToImage(RC(:,2) - RCtr(:,2), RC(:,1) - RCtr(:,1));
dcolor = reshape(dcolor, [n n 3]);