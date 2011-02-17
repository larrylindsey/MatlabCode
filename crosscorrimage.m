function imax = crosscorrimage(imlt, imrt, lim)
if nargin < 3
    lim = 500;
end

step = 10;
iter = step:step:lim;

out = zeros(1, length(iter));

for i_iter = 1:length(iter)
    ii = iter(i_iter);
    ltend = imlt(1:2:end,(end-ii+1):end);
    rtend = imrt(1:2:end,1:ii);
    ltend = ltend(:) - mean(ltend(:));
    rtend = rtend(:) - mean(rtend(:));
    out(i_iter) = sum(ltend(:) .* rtend(:)) / ii;
end

[corr, imax] = max(out);
imax = iter(imax);
newiter = (-step:1:step) + imax;
out2 = zeros(1,1+2*step);
for i_iter = 1:length(newiter)
    ii = newiter(i_iter);
    ltend = imlt(:,(end-ii+1):end);
    rtend = imrt(:,1:ii);
    ltend = ltend(:) - mean(ltend(:));
    rtend = rtend(:) - mean(rtend(:));
    out2(i_iter) = sum(ltend(:) .* rtend(:)) / ii;
end
[corr2, imax2] = max(out2);
imax = newiter(imax2);




