function [Q,D]=qscore3(cat,wcat,verbose)
% QSCORE3 3D segmentation quality assesment
% [Q,D]=qscore3(wTEST,wGT)
%   wTEST       - TEST segmentation
%   wGT         - GT segmentation
%
%   Q           - overall difference score in clicks;
%                 assuming 3.6K clicks per hour gives
%                 approximate proofreading time
%   D           - detailed breakdown of difference score
%
%  BY Y. Mishchenko  FEB 2008
%
% See code for additional parameters affecting results of
% filtering insignificant fragments.

if(nargin<3)
  verbose = true;
end;

%PARAMETERS THAT AFFECT TEST->GT CORRESPONDENCE
%insignificance relative represented volume limit for 2D/3D form
athr=0.9;
%insignificance size limit for abs area form
bthr=100;
%use insignificance criterion in 2D form?
insignificant_2D=1;
%use insignificance criterion in 3D form?
insignificant_3D=0;
%use insignificance criterion in abs area form?
insignificant_abs=1;
%force 8-connected topology
% [if off, make sure segmentations are 8-connected]
force8=0;



%CONSTRUCT ALLEGIANCE OF 2D SEGMENTS
if(verbose)
  fprintf('Preparing to evaluate...\n');
end;
%have TEST/GT segmentation conform with 8-connected neighb
if(force8)
  for i=1:size(cat,3) cat(:,:,i)=imerode(cat(:,:,i),ones(2)); end
  for i=1:size(wcat,3) wcat(:,:,i)=imerode(wcat(:,:,i),ones(2)); end
end

conn=zeros(3,3,3);
conn(:,:,2)=[0 1 0; 1 1 1; 0 1 0];
%slc is 2D segmentation in TEST
slc=uint32(bwlabeln(cat>0,conn));
%slcGT is 2D segmentation in GT
slcGT=uint32(bwlabeln(wcat>0,conn));

%PARSE 2D SEGMENTS: ALLEGIANCE/SIZE/STATS
ss=regionprops(slc,'PixelIdxList','Area');
for i=1:length(ss)
  %2D segments Z-position
  idx=ss(i).PixelIdxList(1);
  [x,y,z]=ind2sub(size(wcat),idx);
  ss(i).Z=z;

  %this obtains regions in GT that cross 2D segment #i
  tmp1=double(sort(wcat(ss(i).PixelIdxList(:))));
  tmp2=diff([tmp1;max(tmp1)+1]);
  ssum=diff(find([1;tmp2]));
  sidx=tmp1(tmp2>0);
  %areas
  ssumGT=ssum(sidx>0);
  %list of ids
  ccidxGT=uint32(sidx(sidx>0));

  %this obtains regions in GT 2D that cross 2D segment #i
  tmp1=double(sort(slcGT(ss(i).PixelIdxList(:))));
  tmp2=diff([tmp1;max(tmp1)+1]);
  ssum=diff(find([1;tmp2]));
  sidx=tmp1(tmp2>0);
  %areas
  ssumGT2D=ssum(sidx>0);
  %list of ids
  ccidxGT2D=uint32(sidx(sidx>0));

  %this obtains regions in TEST that cross 2D segment #i
  tmp1=double(sort(cat(ss(i).PixelIdxList(:))));
  tmp2=diff([tmp1;max(tmp1)+1]);
  ssum=diff(find([1;tmp2]));
  sidx=tmp1(tmp2>0);
  %areas
  ssumTEST=ssum(sidx>0);
  %list of ids
  ccidxTEST=uint32(sidx(sidx>0));

  %we don't evaluate drawing ==> each 2D segment
  % is forced into 1-1 relationship with all
  % GT, GT 2D and TEST
  if(isempty(ssumGT)) ss(i).GT=0; ss(i).AGT=0;
  else
    [a,b]=max(ssumGT);
    ss(i).GT=ccidxGT(b);
    ss(i).AGT=a;
  end
  if(isempty(ssumGT2D)) ss(i).GT2D=0; ss(i).AGT2D=0;
  else
    [a,b]=max(ssumGT2D);
    ss(i).GT2D=ccidxGT2D(b);
    ss(i).AGT2D=a;
  end
  if(isempty(ssumTEST)) ss(i).TEST=0; ss(i).ATEST=0;
  else
    [a,b]=max(ssumTEST);
    ss(i).TEST=ccidxTEST(b);
    ss(i).ATEST=a;
  end
end

if(insignificant_2D)
  if(verbose)
    fprintf('Determining 2D significance...');
  end;
  %determine significance of 3D regions using 90% coverage in EACH slice
  ssGT2D=regionprops([ss.GT2D],'Area','PixelIdxList');
  unqGT2D=find([ssGT2D.Area]>0);

  ids=[];
  for i=unqGT2D
    ccidx=ssGT2D(i).PixelIdxList;
    ssum=[ss(ccidx).AGT2D];
    if(isempty(ssum)) continue; end

    %determine the smallest number required
    % to provide 90% coverage for 2D GT segment
    [ssum,sidx]=sort(ssum,'descend'); ccidx=ccidx(sidx);
    ssum=cumsum(ssum);
    idx=find(ssum<athr*ssum(end));
    if(isempty(idx)) idx=1; end

    %determine which TEST labels are involved
    % in furnishing representation for 2D GT #i
    ccidx=ccidx(idx);
    ids=union(ids,[ss(ccidx).TEST]);
  end

  idx=1:length(ss);
  idsTEST=[ss(idx).TEST];
  idx=idx(~ismember(idsTEST,ids));
  for j=idx ss(j).Area=0; end
  if(verbose)
    fprintf('discarded %i.\n',length(idx));
  end;
end
if(insignificant_3D)
  %determine significance of 2D segments using 90% coverage in 3D
  if(verbose)
    fprintf('Determining 3D significance...');
  end;
  ssGT=regionprops([ss.GT],'Area','PixelIdxList');
  unqGT=find([ssGT.Area]);

  idsTEST=[ss.TEST];
  idsAREA=[ss.AGT];
  ss1=regionprops(idsTEST,'Area','PixelIdxList');
  ccidx=unique(idsTEST(:))';
  for i=ccidx
    ss1(i).cumArea=sum(idsAREA(ss1(i).PixelIdxList));
  end

  ids=[];
  for i=unqGT
    idx=ssGT(i).PixelIdxList;

    ccidx=unique([ss(idx).TEST]);
    ssum=[ss1(ccidx).cumArea];
    if(isempty(ssum)) continue; end

    %determine the smallest number required
    % to provide 90% coverage for 3D #i
    [ssum,sidx]=sort(ssum,'descend'); ccidx=ccidx(sidx);
    ssum=cumsum(ssum);
    ind=find(ssum<athr*ssum(end));
    if(isempty(ind)) ind=1; end

    %these are insignificant 3D regions
    ccidx=ccidx(ind)';
    ids=union(ids,ccidx);
  end
  idx=1:length(ss);
  idsTEST=[ss(idx).TEST];
  idx=idx(~ismember(idsTEST,ids));
  for j=idx ss(j).Area=0; end
  if(verbose)
    fprintf('discarded %i.\n',length(idx));
  end
end
if(insignificant_abs)
  %DISCARD SMALLER-SIZE 2D FRAGMENTS
  if(verbose)
    fprintf('Determining abs significance...');
  end
  idx=1:length(ss);
  idx=idx([ss.ATEST]<bthr);
  for j=idx ss(j).Area=0; end
  if(verbose)
    fprintf('discarded %i.\n',length(idx));
  end
end

%cleanup
ss=ss([ss.Area]>0);
clear ssumTEST ssumGT ccidxTEST ccidxGT idsTEST idsAREA ids
clear tmp1 tmp2 ssum sidx ccidx n ind idx ss2 j ind idx


%CONSTRUCT SCORE
%this counts how many regions had bee created in
% C={a CROSS b for all a IN TEST AND b IN GT}
regionsC=0;
%this counts split cost in terms of 2D reassignment ops
splits2D=0;
%this counts split cost in terms of 3D split ops
splits3D=0;
ss1=regionprops([ss.TEST],'Area','PixelIdxList');
for i=find([ss1.Area])
  idx=ss1(i).PixelIdxList;

  idsGT=[ss(idx).GT];
  %this obtains list of GT segments in TEST segment, and their Z-spans
  tmp1=double(sort(idsGT(:)));
  tmp2=diff([tmp1;max(tmp1)+1]);
  numGT=diff(find([1;tmp2]));
  idsGT=tmp1(tmp2>0);
  %Z-spans
  numGT=numGT(idsGT>0);
  %list of ids
  idsGT=uint32(idsGT(idsGT>0));
  if(isempty(numGT)) continue; end
  if(length(numGT)>1)
    i=i;
  end

  %number of GT-TEST subregions created
  regionsC=regionsC+length(idsGT);
  %number of 2D TEST that need to be reassigned
  splits2D=splits2D+sum(numGT)-max(numGT);
  %number of 3D fragments of TEST that need to be reassigned
  splits3D=splits3D+length(idsGT)-1;
end

%this counts links costs in terms of 3D link ops
links3D=regionsC-length(unique([ss.GT]));

Q=links3D+splits2D;

D=[];
D.splits2D=splits2D;
D.links3D=links3D;
D.splits3D=splits3D;
