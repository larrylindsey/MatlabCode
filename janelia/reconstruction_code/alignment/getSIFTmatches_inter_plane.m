function matches = getSIFTmatches_inter_plane(frames, descriptors)
%[frames,descriptors,matches]=getSIFTmatches_inter_plane(frames, descriptors, alpha);
% Construct SIFT landmarks and establish matches across the section/plane
% THIS REQUIRES SIFT PACKAGE BY Andrea Vedaldi MAKE SURE IT IS IN THE PATH
%
% INPUTS:
%       frames  is {planes x tiles} cell-array of 4xN SIFT landmarks
%       descriptors is {planes x tiles} cell-array of 128xN SIFT
%              descriptors (FIY)
%       alpha is adjacency matrix, this is 4D array
%             {plane x tiles, plane x tiles} == 1 where two tiles overlap
% OUTPUTS:
%       matches is {planes x tiles x tiles} 3D-cell-array
%           of 2xM match indices into frames
%
% Composed based on A. Vedaldi package by Y.Mishchenko Chklovskii Lab 2008
% v1  09182008  splitting getSIFTmatches into getSIFTframes and
%                 getSIFTmatches - Shiv N. Vitaladevuni, JFRC, HHMI
%

%INTERNAL SWITCHES
v=0;        %set this to make verbous

% -- MATCH ACCROSS PLANES
fprintf('Calculating matches accross planes...\n');
plane_descriptors=cell(1,size(frames,1));
plane_references=cell(1,size(frames,1));
plane_tiles=cell(1,size(frames,1));
for plane=1:size(frames,1)
  for tile1=1:size(frames,2)
    if(isempty(frames{plane,tile1}))
      continue;
    end

    %collect plane-descriptors for matching accross planes
    plane_descriptors{plane}=[plane_descriptors{plane},descriptors{plane,tile1}];
    plane_references{plane}=[plane_references{plane},1:size(descriptors{plane,tile1},2)];
    plane_tiles{plane}=[plane_tiles{plane},repmat(tile1,1,size(descriptors{plane,tile1},2))];
  end
end

for plane1=1:size(frames,1)-1
  fprintf('%d ', plane1);
  plane2=plane1+1;

  %match all descriptors from plane1 onto all descriptors in plane2
  tic
  if(isempty(plane_descriptors{plane1}) || isempty(plane_descriptors{plane2}))
    continue;
  end
  M=siftmatch(plane_descriptors{plane1},plane_descriptors{plane2});
  if(v==1)
    toc;
  end

  %sort out matches back into relevant tiles
  for tile1=1:size(frames,2)
    if(isempty(frames{plane1,tile1}))
      continue;
    end
    for tile2=1:size(frames,2)
      if(isempty(frames{plane2,tile2}))
        continue;
      end
      %identify segment of M that pertains to these two tiles!!!
      idx=find(plane_tiles{plane1}(M(1,:))==tile1 & plane_tiles{plane2}(M(2,:))==tile2);

      matches{plane1,tile1,plane2,tile2} = ...
        [plane_references{plane1}(M(1,idx));plane_references{plane2}(M(2,idx))];
    end
  end

  %SYMMETRIZE ARRAY OF MATCHES -- note switched order in M([2,1],:)
  M=M([2,1],:);
  %sort out matches back into tiles
  for tile1=1:size(frames,2)
    for tile2=1:size(frames,2)
      %identify segment of M that pertains to these two tiles!!!
      idx=find(plane_tiles{plane2}(M(1,:))==tile1 & plane_tiles{plane1}(M(2,:))==tile2);

      matches{plane2,tile1,plane1,tile2} = ...
        [plane_references{plane2}(M(1,idx));plane_references{plane1}(M(2,idx))];
    end
  end
end
fprintf('\n');

return;
end
