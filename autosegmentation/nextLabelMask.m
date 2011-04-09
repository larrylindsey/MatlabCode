function [mask lid] = nextLabelMask(inMask, label, ...
    maskTrans, labelTrans, cstr)
% function [mask lid] = nextLabelMask(inMask, label, ...
%    maskTrans, labelTrans, cstr)
%
% cstr - get an example by calling nextLabelMask with no arguments

if nargin == 0
    cstr.x = [0 1];
    cstr.y = [0 1];
    cstr.iT = .9;
    cstr.eT = .75;
    cstr.imsz = [4096 4096];
    mask = cstr;
    return;
end

x = cstr.x;
y = cstr.y;
if isempty(cstr.imsz)
    imsz = size(inMask);
else
    imsz = cstr.imsz;
end

label = sizeImage(label, imsz);

if ~isempty(maskTrans)
    inMask = applyTransformImage(inMask, maskTrans, x, y, 'nearest');
end

if ~isempty(labelTrans)
    label = applyTransformImage(label, labelTrans, x, y, 'nearest');
end

label = sizeImage(label, size(inMask));
%inMask = sizeImage(inMask, imsz);

[mask lid] = selectLabelsByMask(label, inMask, cstr.iT, cstr.eT);
