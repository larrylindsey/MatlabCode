function [dmag dcolor] = makeDistortionMap(RC, RCtr, t, sz)



[dcolor dmag] = vectorMapToImage(RC(:,2) - RCtr(:,2), RC(:,1) - RCtr(:,1), ...
    t);
dcolor = reshape(dcolor, [sz(1) sz(2) 3]);
dmag = reshape(dmag, sz);