function names = enumerateContours(secdoc)

names = {};


for i_sec = 1:numel(secdoc)
    transforms = secdoc(i_sec).section.Transform;
    for i_t = 1:numel(transforms)
        tform = transforms(i_t);
        for i_c = 1:numel(tform.Contour)
            catName = tform.Contour.name;
            if ~ischar(catName)
                catName = num2str(catName);
            end
            names = unique({names{:} catName});
        end
    end    
end

end