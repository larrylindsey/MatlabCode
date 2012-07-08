function [T fixorder] = rectifyOrder(T, morder, o)

fixorder = listOrder(size(morder,2), o);

check = dist2(morder, fixorder);

[~, reindex] = min(check, [], 2);
T = T(reindex,:);

 
