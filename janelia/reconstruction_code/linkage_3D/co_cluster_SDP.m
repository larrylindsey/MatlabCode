function P = co_cluster_SDP(F)
F = double(F);

P = sdpvar(size(F,1), size(F,1));

constraints = [P>= 0] +  [diag(P)==1] + [0 <= P(:)];

fprintf('optimizing ...\n');
tic
options = sdpsettings('solver','sdpt3');
% options = sdpsettings('solver','sedumi');
solvesdp(constraints, sum(F(:).*P(:)), options);
fprintf('done.\n');
toc

P = double(P);
return;
end
