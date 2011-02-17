function peaks = findPeaks2(in)

if any(isnan(in))
    in(logical(isnan(in))) = min(in);
end

if ndims(in) ~= 2
    error('only works in two dimensions');
end

ind1 = peaks1(in, 1);
ind2 = peaks1(in, 2);

peaks = and(ind1, ind2);

end

function pk = peaks1(in, dim)

if dim == 2
    pk = peaks1(in', 1)';
    return;
end

if any(isnan(in))
    in(logical(isnan(in))) = min(in);
end

inD = diff(in, 1, dim);

signum = sign(inD);
signum(logical(signum == 0)) = 1;

ff = [0; -1; 1];

peakIndicator = conv2(signum, ff, 'same');

pk = logical(peakIndicator == 2);

pk = cat(1, pk, false(1, size(pk, 2)));

end