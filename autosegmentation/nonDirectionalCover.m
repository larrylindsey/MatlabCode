function [pOvlp nOvlp] = nonDirectionalCover(mask1, mask2)

n1 = sum(mask1(:));
n2 = sum(mask2(:));
maskOvlp = mask1 .* mask2;

nOvlp = sum(maskOvlp(:));

if n1 == 0 || n2 == 0
    pOvlp = 0;
else
    pOvlp = nOvlp / sqrt(n1 .* n2);
end

end