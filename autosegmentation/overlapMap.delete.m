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
        
        ns_count = numel(ns_seg);
        
        if ns_count > 0
            
            jmap = zeros(ns_count , l);
            
            jmap(:,1) = jseg;
            jmap(:,2) = ns_seg;
            
            swsel = jmap(:,2) < jmap(:,1);
            jmap(swsel,1:2) = jmap(swsel,[2 1]);
            
            nsegmasks = repmat(nextSlice, [1 1 ns_count ]) == ...
                repmat(reshape(ns_seg, [1 1 ns_count]), [size(nextSlice), 1]);
            
            jmap(:,3:end) = areaMetricsAll(...
                repmat(segmask, [1 1 ns_count]), nsegmasks);
            
%             
%                     for kk = 1:numel(ns_seg)
%                         kseg = ns_seg(kk);
%             
%                         smap = zeros(1, l);
%             
%                         if kseg > jseg
%                             smap(1:2) = [jseg kseg];
%                         else
%                             smap(1:2) = [kseg jseg];
%                         end
%             
%                         nsegmask = nextSlice == kseg;
%             
%                         smap(3:end) = areaMetrics(segmask, nsegmask);
%             
%                         jmap(kk,:) = smap;
%                     end
%             
            map = cat(1, map, jmap);
        end
        
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


function measurements = areaMetricsAll(mask1, mask2)
union = or(mask1, mask2);
intersect = and(mask1, mask2);

s_intersect = doubleSum(intersect);

a = s_intersect ./ doubleSum(union);
b1 = s_intersect ./ doubleSum(mask1);
b2 = s_intersect ./ doubleSum(mask2);
b = max(b1, b2);
c = s_intersect;

measurements = [a b c];
end

function s = doubleSum(mask3d)
s = squeeze(sum(sum(mask3d,1),2));
end
