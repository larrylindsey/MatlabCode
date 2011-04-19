function [e v D] = doNCuts(W)
% [e v D] = doNCuts(W)
% Performs an ncuts calculation on the graph-connection weight matrix W
% 
% For instance, to select the primary segmentation, do
% D * v > 0
% Use the column corresponding to the smallest nonzero eigenvalue to select
% the segmentation groups.

n = 4;

if size(W,1) ~= size(W,2)
    error('W should be square');
end

s = size(W, 1);


D = sparse(1:s, 1:s, sum(W,1), s, s);
Dinv = sparse(1:s, 1:s, 1./sum(W,1), s, s);


M = sqrt(Dinv) * (D - W) * sqrt(Dinv);

[v d] = eigs(M, n, 'SM');
e = d(logical(eye(n)));

