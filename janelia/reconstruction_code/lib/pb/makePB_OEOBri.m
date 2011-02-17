function [pb,pbsoft] = makePB_OEOBri( modelBest, img, rFitParab )
% OE + Oriented Brightness
if nargin<3, rFitParab = 2.5; end

pbRaw = applyPBmodelOEOBri( modelBest, img, rFitParab );
pb = max( pbRaw, [], 3 );
if nargout>1
    pbsoft = max( gentleNonmax(pbRaw), [], 3 );
end
