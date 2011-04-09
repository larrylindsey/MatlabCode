function [pOvlp nOvlp] = directionalCover(innerMask, outterMask)

nFrom = sum(innerMask(:));
maskOvlp = innerMask .* outterMask;

nOvlp = sum(maskOvlp(:));
pOvlp = nOvlp / nFrom;

end