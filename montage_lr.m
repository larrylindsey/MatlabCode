function [ltend rtend himax vimax hcorr vcorr] = montage_lr(imlt, imrt, lim)
if nargin < 3
    lim = 500;
end

%Horizontal Montage

%Coarse - first pass
step = 10;
iter = step:step:lim;
out = zeros(1, length(iter));

for i_iter = 1:length(iter)
    ii = iter(i_iter);
    ltend = imlt(2:2:end,(end-ii+1):end);
    rtend = imrt(2:2:end,1:ii);
    ltend = ltend(:) - mean(ltend(:));
    rtend = rtend(:) - mean(rtend(:));
    out(i_iter) = sum(ltend(:) .* rtend(:)) / numel(ltend) / numel(rtend);
end

[vcorr, imax] = max(out);
imax = iter(imax);

%Fine - second pass
iter = (-step:1:step) + imax;
out = zeros(1,1+2*step);
for i_iter = 1:length(iter)
    ii = iter(i_iter);
    ltend = imlt(:,(end-ii+1):end);
    rtend = imrt(:,1:ii);
    ltend = ltend(:) - mean(ltend(:));
    rtend = rtend(:) - mean(rtend(:));
    out(i_iter) = sum(ltend(:) .* rtend(:)) / numel(ltend) / numel(rtend);
end
[vcorr, imax] = max(out);
himax = iter(imax);

%Vertical
%Coarse
chunk_lt = imlt(:,(end-himax+1):end);
chunk_rt = imrt(:,1:himax);
limv = floor(lim/2);
iter = -limv:step:limv;
out = zeros(1, length(iter));

for i_iter = 1:length(iter)
    ii = iter(i_iter);
    ltend = chunk_lt(max(1, 1+ii):min(end, end + ii),2:2:end);
    rtend = chunk_rt(max(1, 1-ii):min(end, end - ii),2:2:end);
    %keyboard
    ltend = ltend(:) - mean(ltend(:));
    rtend = rtend(:) - mean(rtend(:));
    out(i_iter) = sum(ltend(:) .* rtend(:)) / numel(ltend) / numel(rtend);
end
[hcorr, imax] = max(out);
imax = iter(imax);

%Fine
iter = (-step:1:step) + imax;
out = zeros(1,1+2*step);
for i_iter = 1:length(iter)
    ii = iter(i_iter);
    ltend = chunk_lt(max(1, 1+ii):min(end, end + ii),:);
    rtend = chunk_rt(max(1, 1-ii):min(end, end - ii),:);
    %keyboard
    ltend = ltend(:) - mean(ltend(:));
    rtend = rtend(:) - mean(rtend(:));
    out(i_iter) = sum(ltend(:) .* rtend(:)) / numel(ltend) / numel(rtend);
end
[vcorr, imax] = max(out);
vimax = iter(imax);
ii = vimax;
ltend = chunk_lt(max(1, 1+ii):min(end, end + ii),:);
rtend = chunk_rt(max(1, 1-ii):min(end, end - ii),:);
%keyboard



