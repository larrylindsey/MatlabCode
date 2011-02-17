function [n fp] = bestNumPosts(pk, maxN)

npk = numel(pk);

if npk >= maxN    
    fp = linspace(pk(1), pk(end), maxN);
    n = maxN;
    return;
end

e = nan(1, maxN);
pk = sort(pk);

for i_n = npk:maxN
    fp = linspace(pk(1), pk(end), i_n);
    e(i_n) = calcError(fp, pk);
end

n = findPeaks(-e, 1, 'first');
fp = linspace(pk(1), pk(end), n);


end

function eout = calcError(fp, pk)

if size(fp, 2) > size(fp, 1)
    fp = fp';
end

if size(pk, 2) > size(pk, 1)
    pk = pk';
end

 dd = dist2(fp, pk);
 
 [dmin ifpmin] = min(dd, [], 1);
 assign_fp = zeros(size(pk));
 
 [junk ipksort] = sort(dmin);
 
 for i_ps = 1:numel(ipksort)
     i_p = ipksort(i_ps);
     assign_fp(i_p) = ifpmin(i_p);
     dd(assign_fp(i_p), :) = nan;
     [dmin ifpmin] = min(dd, [], 1);
 end

 fp = fp(assign_fp);
 
 eout = rms(fp - pk);
 
end