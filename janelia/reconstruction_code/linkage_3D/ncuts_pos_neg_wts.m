function [V, L, C] = ncuts_pos_neg_wts(W, delta, k, is_verbose)

if(nargin<=2)
  is_verbose = false;
end

A = zeros(size(W));
A(W>0) = W(W>0);
A = A + delta;

R = zeros(size(W));
R(W<0) = -W(W<0);
R = R + delta;

D_A = diag(sum(A,2));
D_R = diag(sum(R,2));

W_eq = A - R + D_R;
D_eq = D_A + D_R;
D_eq = max(D_eq, eye(size(D_eq)));

[V, L] = eigs(W_eq, D_eq, k, 'LA');
if(is_verbose)
  fprintf('W_eq:\n');
  display(W_eq);
  fprintf('D_eq:\n');
  display(D_eq);
  fprintf('V:\n');
  display(V);
  fprintf('L:\n');
  display(L);
end

% x = V.*repmat(abs(diag(L))', [size(A,1), 1]);
x = V;
d = pdist(x, 'euclidean');
C = squareform(d);
if(is_verbose)
  fprintf('C:\n');
  display(C);
end

return;
end
