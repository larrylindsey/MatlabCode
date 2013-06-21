function map = overlapMap(bwl3d)
% map = overlapMap(bwl3d)
%  Given a 3d labeling, creates an overlap map, in [m l]. The first two columns
%   represent label indices, the remaining columns are filled with weight
%   calculations, like overlap area, intersection area, etc.

n = size(bwl3d, 3);
l = 5;

map = [];

for ii = 1:(n - 1)
    slice = bwl3d(:,:,ii);
    nextSlice = bwl3d(:,:,ii+1);
    
    s_seg = unique(slice(:));
    s_seg(s_seg == 0) = [];
    
    for jj = 1:numel(s_seg)
        jseg = s_seg(jj);        
        segmask = slice == jseg;
        
        ns_seg = unique(nextSlice(segmask));
        ns_seg(ns_seg == 0) = [];
        
        jmap = zeros(numel(ns_seg), l);
        
        for kk = 1:numel(ns_seg)
            kseg = ns_seg(kk);
            
            smap = zeros(1, l);
            
            if kseg > jseg
                smap(1:2) = [jseg kseg];
            else
                smap(1:2) = [kseg jseg];
            end
            
            nsegmask = nextSlice == kseg;
            
            smap(3:end) = areaMetrics(segmask, nsegmask);
            
            jmap(kk,:) = smap;
        end
        
        map = cat(1, map, jmap);
        
    end
end

end

function measurements = areaMetrics(mask1, mask2)
union = or(mask1, mask2);
intersect = and(mask1, mask2);

s_intersect = sum(intersect(:));

a = s_intersect / sum(union(:));
b1 = s_intersect / sum(mask1(:));
b2 = s_intersect / sum(mask2(:));
b = max(b1, b2);
c = s_intersect;

measurements = [a b c];
end
