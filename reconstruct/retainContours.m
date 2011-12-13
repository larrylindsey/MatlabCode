function secdoc = retainContours(secdoc, cnames)

parfor idoc = 1:numel(secdoc)
    trans = secdoc(idoc).section.Transform;
    iTransImage = secdoc(idoc).section.transImageIndex;
    for itrans = 1:numel(trans)
        if itrans ~= iTransImage            
            contours = trans(itrans).Contour;
            sel = true(size(contours));
            for icontour = 1:numel(contours)
                cname = contours(icontour).name;
                sel(icontour) =  ~isempty(strsearch(cnames, ...
                    ['^' cname '$']));
            end
            trans(itrans).Contour = contours(sel);            
        end
    end
    secdoc(idoc).section.Transform = trans;
end

end