function peaks = findPeaks(in, varargin)

if any(isnan(in))
    in(logical(isnan(in))) = min(in);
end

inD = diff(in);

signum = sign(inD);
signum(logical(signum == 0)) = 1;

peakIndicator = conv2(signum, [-1 1], 'same');

peaks = find(peakIndicator == 2, varargin{:}) + 1;