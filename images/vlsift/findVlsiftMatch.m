function imatch = findVlsiftMatch(d1, d2, n, type)

if nargin < 4
    mfun = @eucMatch;
    doNorm = false;
else
    if strcmpi(type(1:3), 'euc')
        mfun = @eucMatch;
        doNorm = false;
    else
        mfun = @dotMatch;
        doNorm = true;
    end
end
    

d1 = d1';
d2 = d2';

if doNorm
    d1n = repmat(sqrt(sum(d1.^2, 2)), [1 size(d1, 2)]);
    d2n = repmat(sqrt(sum(d2.^2, 2)), [1 size(d2, 2)]);
    d1 = d1 ./ d1n;
    d2 = d2 ./ d2n;
end

imatch = zeros(size(d1, 1), 1);

istart = 1:n:numel(imatch);
istart = istart(1:(end - 1));
ifinish = [istart(2:end) - 1 numel(imatch)];

for ii = 1:numel(istart)
    currsel = istart(ii):ifinish(ii);
    imatch(currsel) = mfun(d1(currsel,:), d2);
end

end

function imatch = eucMatch(d1, d2)
dist = dist2(d1, d2);
[junk imatch] = min(dist, [], 2);%#ok
end

function imatch = dotMatch(d1, d2)
aff = d1 * d2';
[junk imatch] = max(aff, [], 2);%#ok
end