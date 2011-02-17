function peaks = findPeaks3(in)

if any(isnan(in))
    in(logical(isnan(in))) = min(in);
end

if ndims(in) ~= 3
    error('only works in three dimensions');
end

ind1 = peaks1(in, 1);
ind2 = peaks1(in, 2);
ind3 = peaks1(in, 3);

peaks = logical(ind1 .* ind2 .* ind3);

end

function pk = peaks1(in, dim)

if dim ~= 1
    perm = 1:3;
    perm(dim) = 1;
    perm(1) = dim;
    in = permute(in, perm);
    
    pk = peaks1(in, 1);
    
    pk = ipermute(pk, perm);
    return;
end

if any(isnan(in))
    in(logical(isnan(in))) = min(in);
end

inD = diff(in, 1, 1);

signum = sign(inD);
signum(logical(signum == 0)) = 1;

ff = [0; -1; 1];

peakIndicator = convn(signum, ff, 'same');

pk = logical(peakIndicator == 2);

pk = cat(1, pk, false(1, size(pk, 2), size(pk, 3)));

end