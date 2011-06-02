function [e v D] = doNCuts(W, n)
% [e v D] = doNCuts(W, n)
% Performs an ncuts calculation on the graph-connection weight matrix W
% 
% For instance, to select the primary segmentation, do
% D * v > 0
% Use the column corresponding to the smallest nonzero eigenvalue to select
% the segmentation groups.
%
% Returns the first n smalled eigenvalues, or 4 if n is unspecified.

if nargin < 2
    n = 4;
end

if size(W,1) ~= size(W,2)
    error('W should be square');
end

s = size(W, 1);


D = sparse(1:s, 1:s, sum(W,1), s, s);
Dinv = sparse(1:s, 1:s, 1./sum(W,1), s, s);


M = sqrt(Dinv) * (D - W) * sqrt(Dinv);

M = enforceSymmetric(M);

[v d] = eigs(M, n, 'SM');
e = d(logical(eye(n)));

