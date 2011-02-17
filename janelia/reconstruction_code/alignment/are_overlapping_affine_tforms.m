function are_overlapping = are_overlapping_affine_tforms(...
  tform_1, tform_2, size_1, size_2)
% finds whether two images with affine transforms tform_1 and tform_2 are
% overlapping.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

global config_global

n_points = 11;
zs = zeros(1, n_points);
os = ones(1, n_points);
ps = 0:1/(n_points-1):1;
psr = flipdim(ps,2);
%Get bounding rectangles for tiles
X1 = [size_1(1) * [zs, ps, os, psr]; ...
  size_1(2) * [ps, os, psr, zs]];
X2 = [size_2(1) * [zs, ps, os, psr]; ...
  size_2(2) * [ps, os, psr, zs]];

X1_t = tform_1* [X1; ones(1, size(X1,2))];
X2_t = tform_2* [X2; ones(1, size(X2,2))];

if(isfield(config_global, 'debug_verbose_figures') && ...
    config_global.debug_verbose_figures)
  figure(456);
  plot(X1_t(1,:), X1_t(2,:), '.');
  hold on;
  plot(X2_t(1,:), X2_t(2,:), 'r.');
  hold off;
  title('are_overlapping_affine_tforms');
end

are_overlapping = ...
  max(inpolygon(X1_t(1,:), X1_t(2,:), X2_t(1,:), X2_t(2,:)))>0;

return;
end
