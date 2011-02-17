function P = co_cluster_SDP_cvx(F)
F = double(F);

n_node = size(F,1);

cvx_begin SDP
  variable P(n_node, n_node) symmetric
  variable B(n_node, n_node) symmetric
  minimize( sum(sum(F.*P)) )
  subject to
    P == semidefinite(n_node);
    diag(P) == ones(n_node,1);
    P(:) >= 0;
    P(:) <= 1;
cvx_end

P = double(P);
return;
end
