function gridMatchExampleFigures(n, rdist, sh)

t = linspace(-1,1,n);
XY = gridRC(t,t);
r = repmat(sqrt(sum(XY.^2, 2)), [1 2]);

tr = identityTransform(@taylorMat, 1);
tr.T(2) = sh;
XY = applyTransform(XY, tr);

XYtr = XY .* (1 + rdist * r);

XYaff = affineAlign(XY, XYtr);

clf;
hold on;
%plot(XYaff(:,1), XYaff(:,2), 'kx', 'LineWidth', 2);
XYaff = reshape(XYaff, [n n 2]);

for ii = 1:n
    plot(squeeze(XYaff(:,ii,1)), squeeze(XYaff(:,ii,2)), 'k', 'LineWidth', 1);
    plot(squeeze(XYaff(ii,:,1)), squeeze(XYaff(ii,:,2)), 'k', 'LineWidth', 1);
end


plot(XYtr(:,1), XYtr(:,2), 'bo', 'LineWidth', 2);

a = max(abs(XYtr(:)));
b = max(abs(XYaff(:)));
a = max(a, b) * 1.1;

axis image;
axis([-1 1 -1 1] * a);

set(gca, 'XTickLabel', '', 'YTickLabel', '')
