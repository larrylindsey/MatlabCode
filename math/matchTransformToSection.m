function tr = matchTransformToSection(tr, imagefiles, secdoc)

if numel(tr) ~= numel(imagefiles)
    error('tr and imageFiles must contain the same number of elements');
end

[secfiles, secindex] = enumerateImages(secdoc);

for ii = 1:numel(tr)
    [~, si] = strsearch(secfiles, imagefiles{ii});
    tr(ii).secIndex = secindex(si);    
end


end
