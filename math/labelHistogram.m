function [h o] = labelHistogram(l, doOffset)

if nargin > 1 && doOffset
    o = min(l(:)) - 1;
    l = l - o;
else
    o = 0;
end

n = max(l(:));

h = zeros(n, 1);

parfor il = 1:n
    h(il) = sum(l(:) == il); %#ok<PFBNS>
end
