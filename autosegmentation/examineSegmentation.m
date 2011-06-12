function examineSegmentation(secdoc, segstack, indices)

sections = [secdoc.section];
secIndices = [sections.index];

ax = zeros(1, 4);

for ii = 1:4
    ax(ii) = subplot(2, 2, ii);
end

allseg = sum(segstack, 3) > 0;
[rr cc] = find(allseg);
minR = min(rr); maxR = max(rr); minC = min(cc); maxC = max(cc);
rangeR = maxR - minR; rangeC = maxC - minC;

minR = minR - rangeR;
maxR = maxR + rangeR;
minC = minC - rangeC;
maxC = maxC + rangeC;
selC = minC:maxC;
selR = minR:maxR;

selC = selC(selC > 0);
selC = selC(selC <= size(allseg, 2));
selR = selR(selR > 0);
selR = selR(selR <= size(allseg, 1));

x = [];
y = [];
% selR = [];
% selC = [];

for ii = 1:numel(indices)
    index = indices(ii);
    sec = secdoc(secIndices == index);
    
    imname = sprintf('seg_cache/image_%04d.png', index);
    fid = fopen(imname, 'r');
    if fid < 0
        if isempty(x)
            [x y] = reconstructDomainBounds(secdoc);
        end
        imindex = sec.section.transImageIndex;
        tr = sec.section.Transform(imindex);
        im = imread(tr.Image.src);
        if size(im, 3) > 1
            im = rgb2gray(im);
        end
        im = im2double(im);
        im = applyTransformImage(im, tr, x, y);
        imwrite(im, imname);
    else
        fclose(fid);
        im = imread(imname);
    end
    
%     oldSelR = selR;
%     oldSelC = selC;
%     [selR selC] = getWindowSel(segstack(:,:,ii));
%     
%     if isempty(selR) || isempty(selC)
%         selC = oldSelC;
%         selR = oldSelR;
%     end
    
    im = im2double(im(selR, selC));
    imc = imread(sprintf('seg_cache/class_%04d.png', index));
    imc = im2double(imc(selR, selC));
    iml = imread(sprintf('seg_cache/label_%04d.png', index));
    iml = iml(selR, selC);
    
    imlrgb = im2double(label2rgb(iml, 'hsv', 'k', 'shuffle'));
    
    imshow(imlrgb .* sizeImage(im, size(imlrgb)), 'Parent', ax(1));
    imshow(.5 * (im + segstack(selR, selC, ii)), 'Parent', ax(2));
    imshow(imc, 'Parent', ax(3));
    imshow(im, 'Parent', ax(4));
    title(sprintf('Index %d', index));
    drawnow
    
    imwrite(.5 * (imlrgb + repmat(segstack(selR, selC, ii), [1 1 3])) .* ...
        sizeImage(im, size(imlrgb)), ...
        sprintf('example_%04d.png', ii));
    
end
end

function [selR selC] = getWindowSel(seg)
[rr cc] = find(seg);
minR = min(rr); maxR = max(rr); minC = min(cc); maxC = max(cc);
rangeR = maxR - minR; rangeC = maxC - minC;

minR = minR - rangeR;
maxR = maxR + rangeR;
minC = minC - rangeC;
maxC = maxC + rangeC;
selC = minC:maxC;
selR = minR:maxR;

selC = selC(selC > 0);
selC = selC(selC <= size(seg, 2));
selR = selR(selR > 0);
selR = selR(selR <= size(seg, 1));
end