function [trinfo] = collectReconstructTransforms(secdoc)

fwfunc = @doCubicTransform;
revfunc = @doInverseCubicTransform;

for i_s = 1:numel(secdoc)    
    sec = secdoc(i_s);
    trans = sec.section.Transform(sec.section.transImageIndex);
    trans.data.u = [0 max(trans.Contour.points(:,1)) * ...
        trans.Image.mag];
    trans.data.v = [0 max(trans.Contour.points(:,2)) * ...
        trans.Image.mag];
    
     trinfo1.tr = maketform('custom', 2, 2, fwfunc, revfunc, trans);
%     trinfo1.u = udata;
%     trinfo1.v = vdata;
     trinfo1.src = trans.Image.src;
     trinfo1.mag = trans.Image.mag;
    
    if i_s == 1
        trinfo = repmat(trinfo1, [1 numel(secdoc)]);
    else
        trinfo(i_s) = trinfo1;
    end
end


end