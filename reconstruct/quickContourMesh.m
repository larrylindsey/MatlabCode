function [tri pts] = quickContourMesh(pts1, pts2, z)

if numel(pts2) > numel(pts1)
    [tri pts] = quickContourMesh(pts2, pts1, z([2 1]));
    return;
end
% 
% hold on;
% plot(pts1(:,1), pts1(:,2), 'b');
% plot(pts2(:,1), pts2(:,2), 'r');


%pts1rad = getAngle(pts1);
%pts2rad = getAngle(pts2);
pts1rad = pts1;
pts2rad = pts2;

dmap = dist2(pts1rad, pts2rad);

[junk i1] = min(dmap, [], 1);
[junk i2] = min(junk);

i1 = i1(i2);

travRecord1 = true(size(pts1, 1), 1);
travRecord2 = true(size(pts2, 1), 1);

preTri = [];

while travRecord1(i1) || travRecord2(i2)
    travRecord1(i1) = false;
    travRecord2(i2) = false;
    
%     plot(pts1(i1,1), pts1(i1,2), 'bx');
%     plot(pts2(i2,1), pts2(i2,2), 'rx');
    
    [iTri, i1, i2] = helperFunction(dmap, i1, i2);

    preTri = cat(2, preTri, iTri);
%     pause;
end

preTri = cat(2, preTri, helperFunction(dmap, i1, i2));

z1 = z(1) * ones(size(pts1, 1), 1);
z2 = z(2) * ones(size(pts2, 1), 1);

pts = cat(1, cat(2, pts1, z1), cat(2, pts2, z2));
tri = preTri(:,:,1);
tri = tri + size(pts1, 1) * (preTri(:,:,2) - 1);

tri = tri';

end

function [iTri i1 i2] = helperFunction(dmap, i1, i2)

i1n = i1 + 1;
i2n = i2 + 1;

if i1n > size(dmap, 1)
    i1n = 1;
end

if i2n > size(dmap, 2)
    i2n = 1;
end

if dmap(i1n, i2) < dmap(i1, i2n)
    iTri = [i1; i2; i1n];
    iTri = cat(3, iTri, [1; 2; 1]);
    i1 = i1n;
else
    iTri = [i1; i2; i2n];
    iTri = cat(3, iTri, [1; 2; 2]);
    i2 = i2n;
end


end

function rad = getAngle(pts)

ctr = mean(pts, 1);

pts = pts - repmat(ctr, [size(pts, 1), 1]);

ptsnorm = repmat(sqrt(sum(pts.^2, 2)), [1 2]);
rad = pts ./ ptsnorm;

%rad = pts;
end