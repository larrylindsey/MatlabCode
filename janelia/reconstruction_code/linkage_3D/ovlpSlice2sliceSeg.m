function ovlpList = ovlpSlice2sliceSeg( seg )
    
% overlaps
ovlpList=struct('sliceLo',{},'loIdx',{},'upIdx',{},'pxls',{});
clear ovlp %will add to ovlpList

if iscell( seg )
    [ny,nx]=size( seg{1} );
    nz=length( seg );
else
    [ny,nx,nz]=size( seg );
end

for is=1:nz-1
    sLo=slice(is);
    sUp=slice(is+1);
    lo=regionprops(sLo,'PixelIdxList');
    up=regionprops(sUp,'PixelIdxList');

    ovlps=[ sLo(:)  sUp(:) ];
    nonborders = ovlps(:,1)>0 & ovlps(:,2)>0;
    ovlps = ovlps( nonborders, :);
    ovlps = unique(ovlps, 'rows');
    
    for io=1:size(ovlps,1)
            ovlp.sliceLo = is;
            ovlp.loIdx = ovlps(io,1);
            ovlp.upIdx = ovlps(io,2);
            ovlp.pxls = intersect( lo(ovlp.loIdx).PixelIdxList, up(ovlp.upIdx).PixelIdxList );
            
            % append to ovlp list
            ovlpList(end+1) = ovlp;
            % back ref ? ovlpIdx( ilo, iup ) = length(ovlpList); % idx of new member
    end
end

%-----
    function s=slice(iz)
        if iscell( seg )
            s=seg{iz};
        else
            s=seg(:,:,iz);
        end
    end

end