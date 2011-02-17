function NImage = NormalizeWithOverlap(Image,PatchSize)
  PImage = padarray(Image,[PatchSize PatchSize],'circular');
  NImage = zeros(size(PImage));
  Count = zeros(size(PImage));
  CentInd = OverlapMosaic(size(PImage),PatchSize,1);
%   TileInd = PatchInd(CentInd,PatchSize,MatrixSize(1));
  for ii=1:length(CentInd)
    TI = PatchInd(CentInd(ii),PatchSize,size(PImage,1));
    NPatch = PImage(TI);
    NPatch = NPatch - mean2(NPatch);
    NPatch = NPatch./std2(NPatch);
    NImage(TI) = NImage(TI) + NPatch;
    Count(TI) = Count(TI) + 1;
  end
  NImage = NImage./Count;
  NImage = NImage((PatchSize+1):(size(PImage,1)-PatchSize),(PatchSize+1):(size(PImage,2)-PatchSize));
end