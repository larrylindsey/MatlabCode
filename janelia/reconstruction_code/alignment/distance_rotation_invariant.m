function D = distance_rotation_invariant(A, B)
% D = distance_rotation_invariant(A, B)
% Assuming B is a rotation of A, find out the errors in pairwise distances
%
% Shiv N. Vitaladevuni
% 

pairwise_diff_A = pdist2(A');
pairwise_diff_B = pdist2(B');

d = abs(pairwise_diff_A - pairwise_diff_B);

D = zeros(size(A,2));
index = 1;
for i = 1:size(A,2)-1
  D(i,i+1:end) = d(index:(index+size(A,2)-i-1));
  index = index + size(A,2)-i;
end
D = max(D,D');

return
end
