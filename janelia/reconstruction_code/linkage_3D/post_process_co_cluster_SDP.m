function PP = post_process_co_cluster_SDP(P)

P = double(P);

Q = chol(P)';

Q_t = zeros(size(Q));
for u = 1:size(Q,1)
  [junk, v] = max(Q(u,:));
  Q_t(u,v) = 1;
end

PP = Q_t*(Q_t');

return;
end
