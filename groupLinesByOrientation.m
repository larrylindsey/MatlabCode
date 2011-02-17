function groups = groupLinesByOrientation(lines)

n = numel(lines);

orientation = zeros(1, n);
tolerance = 8 * 2 * pi / n;
nBins = 360;

for ii = 1:n
    pt1 = lines(ii).point1;
    pt2 = lines(ii).point2;
    xy = pt2 - pt1;
    orientation(ii) = atan2(xy(1), xy(2));    
end

orientationPi = mod(orientation, pi);
orientationHalfPi = mod(orientation, pi / 2);

[hh xx] = hist(orientationHalfPi, nBins);

[junk, iO] = max(hh);
oHalfPi = xx(iO);

totalPopulationIndices = logical(abs(orientationHalfPi - oHalfPi)...
    < tolerance);

totalPopulation = lines(totalPopulationIndices);
tPPiOrientation = orientationPi(totalPopulationIndices);

firstPopulationIndices = logical(abs(tPPiOrientation - oHalfPi) ...
    < tolerance);
secondPopulationIndices = not(firstPopulationIndices);

firstPopulation = totalPopulation(firstPopulationIndices);
secondPopulation = totalPopulation(secondPopulationIndices);

groups = {firstPopulation secondPopulation};

end