function stackOut = autoSegmentReconstruct(secdoc, cname, index, control)

if nargin < 4
    control = defaultControl;
    if nargin == 0
        stackOut = control;
        return;
    else
        control.cstr.imsz = [];
        [x y] = reconstructDomainBounds(secdoc);
        control.cstr.x = x;
        control.cstr.y = y;
        control.cstr.eT = .5;
        control.cstr.iT = .8;
    end
end


if isempty(control.classImage)
    classImFun = @defaultClassImFun;
else
    classImFun = control.classImage;
end

if isempty(control.classThreshold)
    classTh = @defaultClassThreshold;
else
    classTh = control.classThreshold;
end

cstr = control.cstr;


% Get first match mask
annoMask = getRawImageContourSlice(secdoc, index, cname);
classOrig = classImFun(secdoc, index);
labelOrig = fijiTrainableSegmentationAnnotation(classTh(classOrig));
trOrig = secdocTrans(secdoc, index);
maskOrig = nextLabelMask(annoMask, labelOrig, trOrig, trOrig, cstr);

%stackOut = maskOrig;
stackOut = [];
nextMask = maskOrig;

%cstr.imsz = size(maskOrig);

while any(nextMask(:))
    stackOut = cat(3, stackOut, nextMask);
    currMask = nextMask;
    index = index + 1;
    fprintf('Now segmenting index %d\n', index);
    class = classImFun(secdoc, index);
    label = fijiTrainableSegmentationAnnotation(classTh(class));
    tr = secdocTrans(secdoc, index);
    nextMask = nextLabelMask(currMask, label, [], tr, cstr);    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function classIm = defaultClassImFun(section, index)
t = secdocTrans(section, index);
fname = ['class_' t.Image.src];
classIm = im2double(imread(fname));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bw = defaultClassThreshold(imclass)
bw = imclass > .75;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function control = defaultControl
control.classImage = @defaultClassImFun;
control.classThreshold = @defaultClassThreshold;
control.cstr = nextLabelMask;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = secdocTrans(secdoc, index)
if nargin > 1
    section = secdoc([secdoc.index] == index).section;
else
    section = secdoc;
end
imIndex = section.transImageIndex;
t = section.Transform(imIndex);

end