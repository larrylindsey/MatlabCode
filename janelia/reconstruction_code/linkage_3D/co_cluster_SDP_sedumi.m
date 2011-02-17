function P = co_cluster_SDP_sedumi(F)
F = double(F);

m = size(F,1);
n = numel(F);
F = double(F);
c = [zeros(2*n,1); F(:)];
A = sparse([1:n, 1:n, n+1:2*n, n+1:2*n, 2*n+1:2*n+m], ...
  [1:n, 2*n+1:3*n, n+1:2*n, 2*n+1:3*n, 2*n+((1:m)-1)*m+(1:m)], ...
  [-ones(1,n), ones(1,n), ones(1,n), ones(1,n), ones(1,m)], 2*n+m, 3*n);
b = [zeros(n, 1); ones(n, 1); ones(m,1)];
K.l = 2*n;
K.s = size(F,1);

fprintf('optimizing.\n');
tic
[x,y,info] = sedumi(A,b,c,K);
toc
P = reshape(x(2*n+1:end), size(F));

P = double(P);
return;
end
