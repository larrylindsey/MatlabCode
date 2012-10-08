function [imfile it] = reconstructImFile(secdoc, idx)

if nargin < 2
    sec = secdoc(1);
else
    allidx = [secdoc.index];
    sel = find(allidx == idx, 1);
    if isempty(sel)
        error('Section %d not found', idx);
    end
    sec = secdoc(sel);
end

it = sec.section.transImageIndex;

imfile = sec.section.Transform(it).Image.src;
