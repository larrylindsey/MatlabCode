function [pdf ff] = valuesToPDF(hh)

hh = sort(hh(:));

ff = zeros(size(hh));

ff(:) = linspace(0, 1, numel(hh));

[hh ii] = unique(hh);
ff = ff(ii);

pdf = diffinterp1(ff, hh);