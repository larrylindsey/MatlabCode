function AddConstraint(indices, vals, rslt)
%Adds a constraint to sparse normal equation system
%fprintf('Starting AddConstraint\n');
global alignRHS alignLHS
[junk,n] = size(indices);
for i = 1:n
  for j = 1:n
    ii = indices(i);
    jj = indices(j);
    %fprintf('%d %d %f\n', ii, jj, vals(i)*vals(j));
    alignLHS(ii,jj) = alignLHS(ii,jj) + vals(i)*vals(j);
  end
  alignRHS(ii) = alignRHS(ii) + vals(i)*rslt;
end
end