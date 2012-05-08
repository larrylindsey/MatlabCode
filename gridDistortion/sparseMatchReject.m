function [f_a, f_b] = sparseMatchReject(f_a, f_b, k, m)
% function [f_a, f_b] = sparseMatchReject(f_a, f_b, k, [m])
%   Intended for use in matching the location of a small patch in a larger
%   image field, using sift. Rejects spurious matches in the larger field
%   as measured by their k-nearest-neighbor density.
%
%   This function specifically targets frames from the vlfeat matlab
%   toolbox, but may be used over any feature-matching frames.
%
% f_a - frames from the larger field, as calculated from a grayscale image
%       using vl_sift, then matched using vl_ubcmatch. f_a is expected to
%       be a matrix of size [d, n], where d > 2, and n is the number of
%       match points. f_a(1:2,:) must reference the x and y locations of
%       the match point in the image.
% f_b - frames from the smaller field, similar to above, also in [d, n]
% k - the number of nearest neighbors to use for the density calculation.
% m - optionally, the match vector output by vl_ubcmatch. In this case,
%     f_a and f_b may be the frames as output by vl_sift, and this vector is
%     used to select them off.

if nargin > 3
    f_a = f_a(:,m(1,:));
    f_b = f_b(:,m(2,:));
end

dd_a = dist2(f_a(1:2,:)', f_a(1:2,:)');
dd_b = dist2(f_b(1:2,:)', f_b(1:2,:)');

dd_a(logical(eye(size(dd_a)))) = inf;
%Not a typo. dd_a and dd_b should be equal sizes. This will cause an error
%if they aren't.
dd_b(logical(eye(size(dd_a)))) = inf;

[dd_as] = sort(dd_a, 1, 'ascend');
[dd_bs] = sort(dd_b, 1, 'ascend');

pop_density_a = mean(dd_as(1:k,:), 1);
pop_density_b = mean(dd_bs(1:k,:), 1);

sel_a = pop_density_a <= ptile(pop_density_b, .9);

f_a = f_a(:,sel_a);
f_b = f_b(:,sel_a);

