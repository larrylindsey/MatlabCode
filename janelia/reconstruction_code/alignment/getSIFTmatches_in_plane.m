function matches = getSIFTmatches_in_plane(frames, descriptors, alpha)
%[frames,descriptors,matches]=getSIFTmatches_in_plane(frames, descriptors, alpha);
% Construct SIFT landmarks and establish matches within the section/plane
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

if(nargin<3)
  alpha=[];
end

%COMPUTE MATCHES
%MATCHES NEED TO BE COMPUTED DIFFERENTLY WITHIN PLANES AND ACCROSS PLANES
% -- MATCH WITHIN PLANES
fprintf('Calculating matches within planes...\n');
matches=cell([size(frames),size(frames,2)]);
for plane=1:size(frames,1)
  fprintf('%d: ', plane);
  for tile1=1:size(frames,2)
    fprintf('%d ', tile1);
    if(isempty(frames{plane,tile1}))
      continue;
    end

    for tile2=1:size(frames,2)
      if(isempty(frames{plane,tile2}))
        continue;
      end
      if(tile1==tile2)
        continue;
      end
      if(~isempty(alpha) && alpha(plane,tile1,plane,tile2)==0)
        continue;
      end

      %SYMMETRIZE ARRAY OF MATCHES -- note switched order in matches{}([2,1],:)
      if(~isempty(matches{plane,tile2,tile1}))
        matches{plane,tile1,tile2}=matches{plane,tile2,tile1}([2,1],:);
      else
        tic;
        matches{plane,tile1,tile2} = ...
          siftmatch(descriptors{plane,tile1},descriptors{plane,tile2});
      end
    end
  end
  fprintf('\n');
end

return;
end
