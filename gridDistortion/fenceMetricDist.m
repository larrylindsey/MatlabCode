function measure = fenceMetricDist(m_in, samples, L)
% m_in selection
% samples in [n 2] - first peak location, second label
% L [m 2] - first label, second hist value

n = size(samples, 1);
m = size(L, 1);
uselabels = zeros(n, 1);
bwline = false(m, 1);


% Create the boolean vector to be used for the distance transform later.
% We get a 1 wherever there is a peak in the sample set.
for i_n = 1:n
    bwline(samples(i_n, 1)) = true;
    %Figure out which labels we're using
    s = samples(i_n, 1);
    uselabels(i_n) = L(s, 1);    
end

% %Figure out which labels we're using
% for i_m = 1:n    
%     s =samples(i_m, 1);
%     uselabels(i_m) = L(s, 1);    
% end

% Uniquify the used labels.  Is this necessary?
uselabels = unique(uselabels);

% Yeah, probably shouldn't reuse variable names.
s = false(m, 1);

% Create a vector that we can use to select off the histogram values
% corresponding to the labels in use.
for i_m = 1:numel(uselabels)
    sel = logical(L(:, 1) == uselabels(i_m));
    s(sel) = true;
end

% Now, create the selection vector based on the model and the samples
m_select = makeMsel(m_in, s);

bwdist = distanceTransform1d(bwline);
bwdist = bwdist / mean(bwdist);

histo = L(:, 2) .* s;
histo = histo / mean(histo);

%mvect = (bwdist + 1) ./ (histo + 1);

mvect = bwdist;

measure = mean(mvect(m_select));
end

function msel = makeMsel(m, s)

slope = m(1);
offset = m(2);

ii = mod(offset, slope);

minLoc = find(s, 1, 'first');
maxLoc = find(s, 1, 'last');

msel = [];

while ii <= maxLoc
    if ii >= minLoc
        msel = cat(2, msel, ii);
    end
    ii = ii + slope;
end

msel = round(msel);
end