function ms = ptile(m, p, d)
% function ms = ptile(m, p, [d])
%  Finds the given percentile value of a population, along the given
%  dimension
%
%  m - data matrix
%  p - the fractional percentile to find
%  d - the dimension along which to find the percentile. Defaults to the
%      first nonzero dimension if omitted


if nargin < 3
    d = find(size(m) > 1, 1, 'first');
    if isempty(d)
        d = 1;
    end
end

ms = sort(m, d, 'ascend');
pIndex = round(size(m, d) * p);


pvect = 1:numel(size(m));
pvect([1 d]) = [d 1];

ms = permute(ms, pvect);
sz = size(ms);
sz(1) = 1;

ms = ms(pIndex,:);
ms = reshape(ms, sz);
ms = permute(ms, pvect);


