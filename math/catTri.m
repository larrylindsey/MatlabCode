function [tri pts] = catTri(tri1, pts1, tri2, pts2)
tri = cat(1, tri1, tri2 + size(pts1, 1));
pts = cat(1, pts1, pts2);
end