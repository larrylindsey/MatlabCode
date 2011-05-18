function [rmin rmax cmin cmax] = mask2boundingbox(mask)

n = matlabpool('size');
sz = size(mask);

if n == 0
    n = 1;
end

[rmin rmax] = rbound(mask, n, sz);
[cmin cmax] = cbound(mask, n, sz);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rmin rmax] = rbound(mask, n, sz)

rsz = sz(1);

cellRMask = cell(1, n);
rminv = ones(1, n) * sz(1);
rmaxv = ones(1, n);

for p = 1:n
    cellRMask{p} = mask(:,p:n:sz(2));
end

parfor p = 1:n
    for c = 1:size(cellRMask{p},2)
        r = 1;
        while r < rminv(p) && r < rsz
            if cellRMask{p}(r, c)
                rminv(p) = r;
            end
            r = r + 1;
        end
        
        r = rsz;
        while r > rmaxv(p) && r > 1
            if cellRMask{p}(r, c)
                rmaxv(p) = r;
            end
            r = r - 1;
        end
    end
end

rmin = min(rminv);
rmax = max(rmaxv);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cmin cmax] = cbound(mask, n, sz)

csz = sz(2);

cellCMask = cell(1, n);
cminv = ones(1, n) * sz(2);
cmaxv = ones(1, n);

for p = 1:n
    cellCMask{p} = mask(p:n:sz(1),:);
end

parfor p = 1:n
    for r = 1:size(cellCMask{p},1)
        c = 1;
        while c < cminv(p) && c < csz
            if cellCMask{p}(r, c)
                cminv(p) = c;
            end
            c = c + 1;
        end
        
        c = csz;
        while c > cmaxv(p) && c > 1
            if cellCMask{p}(r, c)
                cmaxv(p) = c;
            end
            c = c - 1;
        end
    end
end

cmin = min(cminv);
cmax = max(cmaxv);

end
