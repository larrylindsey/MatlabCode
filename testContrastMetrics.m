function [vim dvim svdim] = testContrastMetrics(im, sq, inc)

sz = size(im);

ith = 1:inc:(sz(1) - sq);
itv = 1:inc:(sz(2) - sq);

vim = zeros(numel(ith), numel(itv));
dvim = vim;
svdim = vim;

sel = (1:sq) - 1;

for ii = 1:numel(ith)
    for jj = 1:numel(itv)
        x = ith(ii);
        y = itv(jj);
        
        imsq = im(sel + x, sel + y);
        
        hdiff = diff(imsq, 1, 1);
        vdiff = diff(imsq, 1, 2);        
        diffim = sqrt(hdiff(:,2:end).^2 + vdiff(2:end,:).^2);        
        
        imsqnorm = imsq / sum(imsq(:));
        
        
        vim(ii, jj) = var(imsq(:));
        dvim(ii, jj) = var(diffim(:));
        svdim(ii,jj) = sum(svd(imsqnorm));
        
    end
end
