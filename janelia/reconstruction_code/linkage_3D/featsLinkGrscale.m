function [logminnorm, lightcos, minavg, litecorr] = featsLinkGrscale(ovlpList,ims)

litecorr=zeros(1,length(ovlpList));
light   =zeros(1,length(ovlpList));
normLo  =zeros(1,length(ovlpList));
normUp  =zeros(1,length(ovlpList));
avgLo   =zeros(1,length(ovlpList));
avgUp   =zeros(1,length(ovlpList));
normmin =zeros(1,length(ovlpList));
avgmin  =zeros(1,length(ovlpList));

is_curr = 1;
imgLo=ims(:,:,is_curr);      
imgUp=ims(:,:,is_curr+1);
for k=1:length(ovlpList), 
    is = ovlpList(k).sliceLo;
    if is~=is_curr
        is_curr=is;
        imgLo=ims(:,:,is_curr);      
        imgUp=ims(:,:,is_curr+1);    
    end

    pxls = ovlpList(k).pxls;
    pxlsLo=imgLo(pxls);
    pxlsUp=imgUp(pxls);
    
    litecorr(k)=corr( pxlsLo, pxlsUp );
    light(k)=sum( pxlsLo .* pxlsUp );
    normLo(k) = norm(pxlsLo);
    normUp(k) = norm(pxlsUp);

    avgLo(k) = sum(pxlsLo) / length(pxls);
    avgUp(k) = sum(pxlsUp) / length(pxls);
    
    minmap=min( imgLo(pxls), imgUp(pxls) );
    normmin(k) = norm(minmap(:));
    avgmin(k) = sum(minmap(:)) / length(pxls);
end

normnorm = normLo .* normUp;
minnorm =min( normLo, normUp );
minavg =min( avgLo, avgUp );
%normnorm( 0==normnorm ) = eps;
lightcos = light ./ normnorm;
lightcos(0==normnorm) = light(0==normnorm) ./ eps;
litecorr(~isfinite(litecorr)) = 0;

logminnorm = log(1+minnorm);
