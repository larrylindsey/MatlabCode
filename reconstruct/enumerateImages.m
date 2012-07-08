function [imnames secindex] = enumerateImages(secdoc)

imnames = cell(numel(secdoc), 1);
secindex = zeros(size(imnames));

for ii = 1:numel(secdoc)
    section = secdoc(ii).section;
    imnames{ii} = section.Transform(section.transImageIndex).Image(1).src;
    secindex(ii) = secdoc(ii).index;
end

end
