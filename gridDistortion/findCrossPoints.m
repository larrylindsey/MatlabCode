function [r c ll] = findCrossPoints(HH, L, method)

if nargin < 3
    centerFunc = @centerOfMass;
else
    switch method
        case 'mass'
            centerFunc = @centerOfMass;
        case 'peak'
            centerFunc = @peak;
        case 'geom'
            centerFunc = @geometricCenter;
        otherwise
            error('Unrecognized method %s', method);
    end
end

n = max(L(:));
r = zeros(1, n);
c = zeros(1, n);
ok = false(1, n);

[rr cc] = meshgrid(1:size(HH,2), 1:size(HH,1));

for ii = 1:n
    labelSelector = logical(L == ii);
    [r(ii) c(ii) ok(ii)] = centerFunc(rr, cc, labelSelector, HH);
end

r = r(ok);
c = c(ok);
ll = 1:n;
ll = ll(ok);


end

function [r c ok] = centerOfMass(rr, cc, labelSelector, HH)

HHMask = HH(labelSelector);
rrMask = rr(labelSelector) .* HHMask;
ccMask = cc(labelSelector) .* HHMask;

% 
% HHMask = HH .* labelSelector;
% rrMask = rr .* HHMask;
% ccMask = cc .* HHMask;

mass = sum(HHMask(:));
rrMass = sum(rrMask(:));
ccMass = sum(ccMask(:));

if mass > 0
    r = rrMass / mass;
    c = ccMass / mass;
    ok = true;
else
    r = 0;
    c = 0;
    ok = false;
end

end

function [r c ok] = peak(rr, cc, labelSelector, HH)
HHMask = HH .* labelSelector;

if max(labelSelector(:)) > 0
    ok = true;
    [i_r i_c] = find(HHMask == max(HHMask(:)));
    r = rr(i_r, i_c);
    c = cc(i_r, i_c);
else
    r = 0;
    c = 0;
    ok = false;
end
end

function [r c ok] = geometricCenter(rr, cc, labelSelector, HH)
if max(labelSelector(:) > 0)
    r = mean(rr(labelSelector));
    c = mean(cc(labelSelector));
    ok = true;
else
    r = 0;
    c = 0;
    ok = false;
end
end