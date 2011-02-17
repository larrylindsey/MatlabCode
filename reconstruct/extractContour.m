function contour = extractContour(secdoc, name)

contour = repmat(struct, 0);
th = [secdoc.section];
th = [th.thickness];
z = cumsum(th) - th(1) + min([secdoc.index]) * median(th);

for i_sec = 1:numel(secdoc)
    transforms = secdoc(i_sec).section.Transform;
    for i_t = 1:numel(transforms)
        tform = transforms(i_t);
        for i_c = 1:numel(tform.Contour)
            if strcmp(name, tform.Contour(i_c).name)
                contour = concatenate(contour, tform.Contour(i_c));
                contour(end).section = secdoc(i_sec).index;
                contour(end).z = z(i_sec);
            end
        end
    end    
end

end

function strout = concatenate(str1, str2)

strout = str1;
if isempty(strout)
    strout = str2;
    return;
end

oldfields = fieldnames(str1);
newfields = fieldnames(str2);

strout(end + 1) = strout(end);

for i_f = 1:numel(newfields)
    name = newfields{i_f};
    strout(end).(name) = str2.(name);
end

nilfields = setdiff(oldfields, newfields);

for i_f = 1:numel(nilfields)
    name = nilfields{i_f};
    strout(end).(name) = [];
end

end
