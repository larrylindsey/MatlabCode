function [pOvlp nOvlp] = nonDirectionalCover(mask1, mask2)

n1 = sum(mask1(:));
n2 = sum(mask2(:));
maskOvlp = mask1 .* mask2;

nOvlp = sum(maskOvlp(:));
pOvlp = nOvlp / sqrt(n1 .* n2);

end