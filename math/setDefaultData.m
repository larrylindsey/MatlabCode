function tr = setDefaultData(tr)
for ii = 1:numel(tr)
    tr(ii).data = struct('n', 32, 'u', [-1 1], 'v', [-1 1]);
end
end
