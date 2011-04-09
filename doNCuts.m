function [e v D] = doNCuts(W)
% [e v D] = doNCuts(W)
% Performs an ncuts calculation on the graph-connection weight matrix W
% 
% For instance, to select the primary segmentation, do
% D * v > 0
% Use the column corresponding to the smallest nonzero eigenvalue to select
% the segmentation groups.


if size(W,1) ~= size(W,2)
    error('W should be square');
end

s = size(W, 1);

D = zeros(size(W));
Dinv = D;

diagsel = logical(eye(s));

D(diagsel) = sum(W, 1);
Dinv(diagsel) = 1./D(diagsel);

M = sqrt(Dinv) * (D - W) * sqrt(Dinv);

[v d] = eig(M);
e = d(diagsel);

