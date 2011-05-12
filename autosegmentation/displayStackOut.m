function displayStackOut(stackout, z)

figure('Renderer', 'OpenGL'); hold on;

for ii = 1:numel(z)
    currstack = stackout(:,:,ii);
    currstack = imclose(currstack, strel('disk', 1));
    P = bwboundaries(currstack);
    for ip = 1:numel(P)
        rc = P{ip};
        
        plot3(rc(:,1), rc(:,2), z(ii) * ones(size(rc(:,1))), 'LineWidth', 2);
    end
end