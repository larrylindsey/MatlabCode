function mout = fenceModelDist(samples, L)
% samples in [n 2] - first peak location, second label
% L [m 2] - first label, second hist value
% ifirst = find(L(:,1) > 0, 1, 'first');
% ilast = find(L(:,1) > 0, 1, 'last');

if isstruct(L)
    L = L.L;
end

%n = size(samples, 1);
% xx = (1:n)';
yy = samples(:, 1);

[n fp] = bestNumPosts(yy, max(L(:,1)) * 2);

mout = [fp(1) fp(2) - fp(1)];


% 
% yy = sort(yy, 'ascend');
% 
% A = cat(2, xx, ones(size(xx)));
% 
% AtA = A' * A;
% [U S V] = svd(AtA);
% Sinv = S;
% S_sel = logical(eye(size(S)));
% Sdiag = S(S_sel);
% Sinv(S_sel) = 1./Sdiag;
% AtAinv = U * Sinv * V';
% 
% mout = AtAinv * A' * yy;
% 




% 
% mout = xx * p(1) + p(2);
% 
% mout = round(mout);
% 
% mout = mout(logical(mout > ifirst));
% mout = mout(logical(mout < ilast));

end