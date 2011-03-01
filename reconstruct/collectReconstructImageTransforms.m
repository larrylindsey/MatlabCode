function imtrans = collectReconstructImageTransforms(secdoc)

imtrans = secdoc(1).section.Transform(secdoc(1).section.transImageIndex);
imtrans = repmat(imtrans, [1 numel(secdoc)]);

for isec = 2:numel(secdoc)
    iti = secdoc(isec).section.transImageIndex;
    currTrans = secdoc(isec).section.Transform(iti);
    ctFNames = fieldnames(currTrans);
    for ifield = 1:numel(ctFNames)
        fn = ctFNames{ifield};
        imtrans(isec).(fn) = currTrans.(fn);
    end
end