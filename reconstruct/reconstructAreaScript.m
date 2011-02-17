function areaplot = reconstructAreaScript(secdoc)

areaplot = zeros(size(secdoc));

for i_sec = 1:numel(secdoc)
    section = secdoc(i_sec).section;
    transform = section.Transform(section.transImageIndex);
    X = transform.Contour.imageDomainTransPoints;
    areaplot(i_sec) = polyarea(X(:,1), X(:,2));
end