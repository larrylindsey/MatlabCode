function [bwl3d, map] = bwlMapEq(bwl3d, eqmap, map)
% [bwl3d, map] = bwlMapEq(bwl3d, eqmap, map)
%    Remaps equivalent segments in bwl3d.
%   bwl3d - 3d label
%   eqmap - like [k n], where k is the number of equivalence relations, and
%           n is at least 2. Columns 3 and above are ignored. The second
%           column must be unique.
%   map - a map as output by overlapMap. The output is reindexed to match
%         the changes to bwl3d.

% Make idx for the sparse matrix
n = max(bwl3d(:));
idx = 1:n;

reqmap = zeros(0,2);
hh = waitbar(0);
% Rectify eqmap
me = max(eqmap(:,2));
for ii = 1:me
    sel = eqmap(:,2) == ii;
    jj = eqmap(sel,1);
    
    nj = numel(jj);
    
    if nj > 1
        eqmap(sel,:) = [];
        remap1 = ones(nj, 1) * jj(1);
        remap2 = jj;
        remap2(1) = ii;
        
        reqmap = cat(1, reqmap, cat(2, remap1, remap2));
        
    end
    waitbar(ii / me, hh);
end

% Make the sparse matrix
Adj = makeNormalizedAdjMat(idx);

% Exponentiate the matrix
AdjM = Adj ^ (2 * size(bwl3d,3)) > 0;


% Run along the rows, zero out the maps for the larger indices.
% If we have
% [1 1 0 1 0
%  1 1 0 1 0
%  0 0 1 0 1
%  1 1 0 1 0
%  0 0 1 0 1 ]
% Then we want
% [1 1 0 1 0
%  0 0 0 0 0
%  0 0 1 0 1
%  0 0 0 0 0
%  0 0 0 0 0]
for ii = 1:size(AdjM,1)
    jj = find(AdjM(ii,:));
    jj(jj == ii) = [];
    AdjM(jj,:) = 0;
    disp(100 * ii / size(AdjM,1));
end

[jj kk] = find(AdjN);

sel = jj ~= kk;
jj = jj(sel);
kk = kk(sel);

idx = 1:n;
for ii = 1:numel(kk)
    idx(jj(ii)) = kk(ii);
end

keyboard

idx = [1 idx + 1];

bwl3d = idx(bwl3d + 1) - 1;

end

function Adj = makeNormalizedAdjMat(idx)
Adj = sparse(idx, 1:numel(idx), ones(size(idx)), numel(idx), numel(idx));
Adj = Adj + Adj';
Adj(Adj > 1) = 1;
nAdj = sum(Adj, 1);
for ii = 1:size(Adj, 2)
    Adj(:,ii) = Adj(:,ii) / nAdj(ii);
end
end



% 
% if nargin < 3
%     map = zeros(1, 2);
% end
% 
% k = size(eqmap, 1);
% n = size(bwl3d, 3);
% 
% [~, isort] = sort(eqmap(:,2), 'descend');
% eqmap = eqmap(isort,1:2);
% eqmapc = {eqmap};
% eqmapc = repmat(eqmapc, [1 k]);
% 
% tsel = zeros(1, n);
% tset = zeros(1, n);
% 
% tbegin = tic;
% 
% parfor jj = 1:n
%     slice = bwl3d(:,:,jj);
%     eqmapj = eqmapc{jj};
%     
% %     bc = 0;
%     
%     for ii = 1:k
%         ts = tic;
%         sel = slice == eqmapj(ii,2);
%         tsel(jj) = tsel(jj) + toc(ts);
%         ts = tic;
%         slice(sel) = eqmapj(ii,1);
%         tset(jj) = tset(jj) + toc(ts);
%         
%         if jj == 1
%             dispstr = sprintf('%d of %d...', ii, k);
% %             fprintf(repmat('\b', [1 bc]));
% %             display(dispstr);
%             system(['echo ' dispstr ' >> eq.log']);
% %             bc = numel(dispstr);
%         end
%     end
%     bwl3d(:,:,jj) = slice;
% end
% 
% toc(tbegin)
% 
% fprintf('selection time: %gs\n', sum(tsel));
% fprintf('set time: %gs\n', sum(tset));
% 
% % 
% % ll = unique(bwl3d(:));
% % ll(ll == 0) = [];
% % ll = sort(ll);
% % 
% % 
% % 
% % for ii = 1:numel(ll)
% %     map2c(map2c == ll(ii)) = ii;
% %     bwl3d(bwl3d == ll(ii)) = ii;
% % end
% % 
% % map(:,1:2) = map2c;
% % 
% % mapsel = map(:,1) < map(:,2);
% % map = map(mapsel,:);
