function [idx c] = hkmeans(x, k, varargin)

c = cell(numel(k), 1);

[idx c0] = kmeans(x, k(1), varargin{:});

c{1} = {c0};

idx = repmat(idx, [1 numel(k)]);
c = repmat(c, [1 numel(k)]);

for isub = 2:numel(k)
    iup = isub - 1;
    lastk = max(idx(:,iup));
    
    ci = cell(lastk, 1);
    
    for jj = 1:lastk
       sel = idx(:,iup) == jj;
       xsel = x(sel,:);
       [idxsub csub] = kmeans(xsel, k(isub), varargin{:});              
       idxsub = idxsub + (jj - 1) * k(isub);     
       
       idx(sel,isub) = idxsub;
       ci{jj} = csub;
    end    
    c{isub} = ci;
end

end
