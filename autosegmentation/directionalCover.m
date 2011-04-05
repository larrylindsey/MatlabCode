function [pOvlp nOvlp] = directionalCover(maskFrom, maskOnto)

nFrom = sum(maskFrom(:));
maskOvlp = maskFrom .* maskOnto;

nOvlp = sum(maskOvlp(:));
pOvlp = nOvlp / nFrom;

end