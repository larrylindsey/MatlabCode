function displaySeriesAlignment(trAlign, imFiles, preTr, delay)

if nargin < 4
    delay = 1;
end

if nargin < 3
    preTr = [];
end

xbound = [inf, -inf];
ybound = xbound;
n = numel(trAlign);

for ii = 1:n
    u = trAlign(ii).data.u;
    v = trAlign(ii).data.v;
    UV = gridRC(u,v);
    XY = applyTransform(UV, trAlign(ii));
    
    xbound(1) = min([xbound(1); XY(:,1)]);
    xbound(2) = max([xbound(2); XY(:,1)]);
    ybound(1) = min([ybound(1); XY(:,2)]);
    ybound(2) = max([ybound(2); XY(:,2)]);
end

lastim = getIm(imFiles{1}, preTr);
lastim = applyTransformImage(lastim, trAlign(1), xbound, ybound);
ih = imshow(repmat(lastim, [1 1 3]));
tic;
for ii = 2:n
    im = getIm(imFiles{ii}, preTr);
    im = applyTransformImage(im, trAlign(ii), xbound, ybound);
    set(ih, 'CData', clipIm(compositeImage(lastim, im)));
    lastim = im;
    title(sprintf('%s-%s', imFiles{ii - 1}, imFiles{ii}));
    drawnow;
    dt = toc;    
    pause(delay - dt);
    tic;
end


end

function im = getIm(name, tr)
MAXSIZE = 512;

im = imread(name);

if size(im,3) > 1
    im = rgb2gray(im);    
end
im = im2single(im);

szim = size(im);

if max(szim) > MAXSIZE
    im = imresize(im, MAXSIZE / max(szim));
end

if ~isempty(tr)
    im = applyTransformImage(im, tr);
end

end

function im = clipIm(im)
im(im > 1) = 1;
im(im < 0) = 0;
end
