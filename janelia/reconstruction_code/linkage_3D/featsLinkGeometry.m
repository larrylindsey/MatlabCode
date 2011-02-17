function [share,width] = featsLinkGeometry(ovlpList,seg)

for k=1:length(ovlpList), 
    ovlp = ovlpList(k);
    areaLo = sum(sum( Slice(ovlp.sliceLo  )==ovlp.loIdx ));
    areaUp = sum(sum( Slice(ovlp.sliceLo+1)==ovlp.upIdx ));
    area=length( ovlp.pxls );
    share(k) = max( area/areaLo, area/areaUp );
    
    if nargout>1     % width
        patch=false(size(Slice(1)));
        patch(ovlp.pxls) = true;
        b=border(patch)&~patch; %outside border
        borderlen=sum(b(:));
        width(k) = area/borderlen;
    end
end

%-----
    function s=Slice(iz)
        if iscell( seg )
            s=seg{iz};
        else
            s=seg(:,:,iz);
        end
    end

end
