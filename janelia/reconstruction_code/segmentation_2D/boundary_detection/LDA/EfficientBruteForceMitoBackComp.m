function Mask = EfficientBruteForceMitoBackComp(Filter,Thresh,Mode)
% function Mask = BruteForceMitoBackComp(Filter,Thresh)
% had to write a back compatible version because cluster does not have
% Matlab 2009 and does not have the function bwconncomp
% Remove mitochondria
  LargeCompSize = 50000;
  SolidityThresh = 0.7;
  SolidityThresh = 0.8;%  Made it more conservative July 7th 2010

  T = Filter<Thresh;
  OpenSize = 2500;
  T = bwareaopen(T,OpenSize);
  se = strel('square',5);
  T = imerode(T,se);
  [S,num] = bwlabel(T);
  Stats = regionprops(S,'Area');
  trid = [0 find([Stats(:).Area]> LargeCompSize)];
  if length(trid)<2
    display('Did not manage to get rid of central component')
    Mask = zeros(size(Filter));
    return
  else
    T = ~ismember(S,trid);
    T = imdilate(T,se);
    T = bwareaopen(T,OpenSize);
    S = bwlabel(T);
    Stats = regionprops(S,'Solidity');
    trid = [0 find([Stats(:).Solidity]<SolidityThresh)];
    T = ~ismember(S,trid);
    se = strel('disk',10);
    T = imclose(T,se);
    T = imdilate(T,ones(3));
    if Mode ~= 3
      Mask = T;
    else
      Mask = false(size(Filter));
    end
    
    if Mode>1
      TT = Filter<Thresh;
      [S,num] = bwlabel(TT);
      Stats = regionprops(S,'Area');
      trid = [0 find([Stats(:).Area]> LargeCompSize)];
      TT = ~ismember(S,trid);
      MinVesicleArea = 20;
      MaxVesicleArea = 600; % Was 400 before
      TT = bwareaopen(TT,MinVesicleArea);
      S = bwlabel(TT);
      Stats = regionprops(S,'Area');
      trid = [0 find([Stats(:).Area]>MaxVesicleArea)];
      TT = ~ismember(S,trid);
      Mask = Mask | TT;
    end
    
  end
end