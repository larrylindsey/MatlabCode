function stats = collectLinearTransformStats(T)

vx = squeeze(T(:,1,:));
vy = squeeze(T(:,2,:));

vv = vx + vy;

innerProduct = dot(vx, vy) ./ (rms(vx, 1) .* rms(vy, 1));
angle = acos(innerProduct);

rotation = atan2(vv(2, :), vv(1,:)) - pi / 4;

mag = T(1,1,:) .* T(2,2,:) - T(1,2,:) .* T(2,1,:);

stretch = rms(vx, 1) ./ rms(vy, 1);

stats = struct('innerProduct', innerProduct(:), 'angle', angle(:),...
    'rotation', rotation(:), 'mag', mag(:), 'stretch', stretch(:));
