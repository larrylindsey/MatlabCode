function [frames,descriptors]=getSIFTframes(I, config)
%[frames,descriptors]=getSIFTframes(I,alpha,config);
% Construct SIFT landmarks
% THIS REQUIRES SIFT PACKAGE BY Andrea Vedaldi MAKE SURE IT IS IN THE PATH
%
% INPUTS:
%       I is {planes x tiles} cell-array of images
%       alpha is adjacency matrix, this is 4D array
%             {plane x tiles, plane x tiles} == 1 where two tiles overlap
%       config is structure of parameters
% OUTPUTS:
%       frames  is {planes x tiles} cell-array of 4xN SIFT landmarks
%       descriptors is {planes x tiles} cell-array of 128xN SIFT
%              descriptors (FIY)
%
% Composed based on A. Vedaldi package by Y.Mishchenko Chklovskii Lab 2008
% v1  09182008  splitting getSIFTmatches into getSIFTframes and
%                 getSIFTmatches - Shiv N. Vitaladevuni, JFRC, HHMI
%

if(nargin<2)
  config=[];
end

if(~isfield(config,'downsampling') || isempty(config.downsampling))
  K=512;
else
  K=config.downsampling;
end

%INTERNAL SWITCHES
v=0;        %set this to make verbose

frames=cell(size(I));       %array of landmarks
descriptors=cell(size(I));  %array of SIFT descriptors
for plane=1:size(I,1)
  fprintf('%d: ', plane);
  for tile=1:size(I,2)
    if(isempty(I{plane, tile}))
      continue;
    end
    fprintf('%d ', tile);
    %DOWNSIZE IMAGES - FACTOR
    S=min(size(I{plane,tile}));
    resampling=round(min(S)/K);

    %DOWNSIZE IMAGES -- sift function takes (!)double(!) image
    H=im2double(imcomplement(imresize(I{plane,tile},1/resampling,'bilinear')));

    %COMPUTE SIFT DESCRIPTORS /ALL ON DEFAULT/
    [frames{plane,tile},descriptors{plane,tile}]=sift(H,'Verbosity',v);
    
    %RESTORE COORDINATES IN FRAMES
    frames{plane,tile}(1:2,:)=frames{plane,tile}(1:2,:)*resampling;
  end
  fprintf('\n');
end

return;
end
