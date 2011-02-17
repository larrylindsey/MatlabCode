function [frames1,matches1,ind1,alpha1,I1]=getSIFTrow(frames,matches,alpha,I)
%[frames,matches,inds,alpha,images]=getSIFTrow(frames,matches,alpha,images)
% Function should be used to collapse 4D representation {planes x tiles,
% planes x tiles} of tiling&registration into 2D vector {tiles',tiles'}
% form that can be used with getSIFTtransform.
%
% INPUTS:
%       frames      is {planes x tiles} cell-array of 4xN SIFT landmarks
%       matches     is {planes x tiles, planes x tiles} cell-array of
%                   2xM match indices into frames
%       alpha       is adjacency 4D array {planes x tiles, planes x tiles}
%       images      is 2D cell-array of images
% OUTPUTS:
%       frames      is 1xK cell-array of 4xN SIFT landmarks
%       matches     is KxK cell-array of 2xM match indices into frames
%       inds        back-reference array of indices
%       alpha       is unwrapped 2D adjacency matrix
%       images      is unwrapped 1xK cell array of images
%
% Composed based on A. Vedaldi package by Y.Mishchenko Chklovskii Lab 2008
if(nargin<3)
  alpha=[];
end

frames1={};
matches1={};
I1={};

%unwrap frames/X
ind1=zeros(size(frames));
cnt1=1;
for plane1=1:size(frames,1)
  for tile1=1:size(frames,2)
    if(isempty(frames{plane1,tile1}))
      continue;
    end
    frames1{cnt1}=frames{plane1,tile1};
    
    if(~isempty(I))
      I1{cnt1}=I{plane1,tile1};
    end
    
    ind1(plane1,tile1)=cnt1;
    cnt1=cnt1+1;    
  end
end


%unwrap matches -- initalize
alpha1=zeros(length(frames1),length(frames1));
for cnt1=1:length(frames1)
  for cnt2=1:length(frames1)
    matches1{cnt1,cnt2}=zeros(2,0);
  end
end

%unwrap matches
for plane1=1:size(frames,1)
  for tile1=1:size(frames,2)
    for plane2=1:size(frames,1)
      for tile2=1:size(frames,2)
        cnt1=ind1(plane1,tile1);
        cnt2=ind1(plane2,tile2);
        if(cnt1==0 || cnt2==0)
          continue;
        end
        matches1{cnt1,cnt2}=matches{plane1,tile1,plane2,tile2};
        if(isempty(matches1{cnt1,cnt2}))
          matches1{cnt1,cnt2}=zeros(2,0);
        end          
        
        if(~isempty(alpha)) 
            alpha1(cnt1,cnt2)=alpha(plane1,tile1,plane2,tile2);
            alpha1(cnt2,cnt1)=alpha1(cnt1,cnt2);
        end
      end
    end
  end
end

if(isempty(alpha))
  alpha1=[];
end
