function names = enumerateContours(secdoc)

p = max(matlabpool('size'), 1);

secdocSplit = cell(1, p);

for ii = 1:p
    secdocSplit{ii} = secdoc(ii:p:end);
end

names = cell(1,p);

parfor ip = 1:p
    names{ip} = {};
    for i_sec = 1:numel(secdocSplit{ip})
        if isfield(secdocSplit{ip}(i_sec).section, 'Transform')
            transforms = secdocSplit{ip}(i_sec).section.Transform;
            for i_t = 1:numel(transforms)
                tform = transforms(i_t);
                for i_c = 1:numel(tform.Contour)
                    catName = tform.Contour(i_c).name;
                    if ~ischar(catName)
                        catName = num2str(catName);
                    end
                    names{ip} = cat(1, names{ip}, catName);
                end
            end
        end
    end
    %names{ip} = unique(names{ip});
end

names = unique(cat(1, names{:}));
end