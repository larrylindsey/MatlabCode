function Sout = enforceSymmetric(Sin)
% function Sout = enforceSymmetric(Sin)
% Copies the upper half matrix into the lower half.

[r c v] = find(Sin);

upperSel = r > c;
diagSel = r == c;

rOut = cat(1, r(upperSel), r(diagSel), c(upperSel));
cOut = cat(1, c(upperSel), c(diagSel), r(upperSel));
vOut = cat(1, v(upperSel), v(diagSel), v(upperSel));

Sout = sparse(rOut, cOut, vOut, size(Sin,1), size(Sin,1));