function [secdoc repSel] = replaceImages(secdoc, oldImages, newImages)

repSel = false(size(oldImages));

for ii = 1:numel(secdoc)
    imIndex = secdoc(ii).section.transImageIndex;
    imTrans = secdoc(ii).section.Transform(imIndex);
        
    images = imTrans.Image;
    
    for jj = 1:numel(images)
        fIndex = strmatch(images(jj).src, oldImages);
        images(jj).src = newImages{fIndex};
        repSel(fIndex) = true;
    end        
    
    secdoc(ii).section.Transform(imIndex).Image = images;
end