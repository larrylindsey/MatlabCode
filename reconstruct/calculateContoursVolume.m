function volumes = calculateContoursVolume(contours)
volumes = zeros(size(contours));
for ii = 1:numel(contours)
    contour = contours{ii};
    slabvols = [contour.slabVol];
    volumes(ii) = sum(slabvols);
end
end
