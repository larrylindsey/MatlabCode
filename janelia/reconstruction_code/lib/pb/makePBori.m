function pbori = makePBori( arrdet, orient, fstd, beta )
% cut out of pbBGTG.m, makePB2
% extended version: cell array of detectors
% each detector, oriented or not, yields one feature
% uses logistic model from segbench
%  Input
%       arrdet -  array of detector responses, y by x by orient or y by x
%       orient  - #orients (or thetalist for back-compat) - common for all oriented detectors
%       fstd, beta - logistic fit
%   Output
%       pb -    map of pb  y*x*theta

[h,w,unused] = size(arrdet{1});
if isscalar(orient), norient=orient;
else norient=length(orient); end %back-compat

pbori = zeros(h,w,norient);
beta = beta ./ [1 fstd];
for i = 1:norient,
  x = ones(h*w,1);
  for idet=1:numel(arrdet)
      if ndims(arrdet{idet})>2
          f = arrdet{idet}(:,:,i);
          x = [x f(:)];
      else
          x = [x arrdet{idet}(:)];
      end
  end
  pbi = 1 ./ (1 + (exp(-x*beta')));
  pbori(:,:,i) = reshape(pbi,[h w]);
end
