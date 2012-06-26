function o = listOrder(k, order)
% function o = listOrder(k, order)
%  Generates a list of polynomial orders in k-space up to the given order
%
%  k - the dimensionality of the space
%  order - the maximal polynomial order requested
%
%  o [l k] - a matrix representing a list of orders per dimension. For
%            example, given two-space and order two, o would be:
%
%     [ 0 0 ;
%       0 1 ;
%       1 0 ;
%       0 2 ;
%       1 1 ;
%       2 0 ]


o = {zeros(1, k)};

for ii = 1:order
    ii1 = ii + 1;
    o{ii1} = [];
    for kk = 1:k
        okk = o{ii};
        okk(:, kk) = okk(:, kk) + 1;
        o{ii1} = cat(1, o{ii1}, okk);
    end
    o{ii1} = unique(o{ii1}, 'rows');
end

o = cat(1, o{:});
