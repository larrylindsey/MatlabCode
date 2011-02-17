function [CentInd,TileInd] = OverlapMosaic(MatrixSize,PatchSize,Overlap)
  if Overlap == 0
    BaseIndx = (PatchSize+1):(2*PatchSize+1):(MatrixSize(1)-(PatchSize));
    BaseIndy = (PatchSize+1):(2*PatchSize+1):(MatrixSize(2)-(PatchSize));
  elseif Overlap == 1
    BaseIndx = (PatchSize+1):PatchSize:(MatrixSize(1)-(PatchSize));
    BaseIndy = (PatchSize+1):PatchSize:(MatrixSize(2)-(PatchSize));
    if(BaseIndx(end) ~= (MatrixSize(1)-PatchSize))
      BaseIndx = [BaseIndx (MatrixSize(1)-PatchSize)];
    end
    if(BaseIndy(end) ~= (MatrixSize(2)-PatchSize))
      BaseIndy = [BaseIndy (MatrixSize(2)-PatchSize)];
    end
% Overcomplicated    
%     Span = MatrixSize(1)-2*PatchSize-1;
%     SpaceVec = PatchSize:-1:1;
%     Vec = Span./SpaceVec;
%     for ii=1:length(Vec)
%       if Vec(ii) == floor(Vec(ii))
%         Space = SpaceVec(ii);
%         break
%       end
%     end
%     BaseInd = (PatchSize+1):Space:(MatrixSize(1)-(PatchSize));
  end
  [I,J] = meshgrid(BaseIndx,BaseIndy);
  CentInd = sub2ind(MatrixSize,I,J);
  CentInd = reshape(CentInd,size(CentInd,1)*size(CentInd,2),1);

  %Potential memory problems
  if nargout == 2
    TileInd = PatchInd(CentInd,PatchSize,MatrixSize(1));
  end
%   TileInd = [];
%    W = (PatchSize+1):(2*PatchSize):(MatrixSize(1)-(PatchSize+1));
%    L = (PatchSize+1):(2*PatchSize):(MatrixSize(2)-(PatchSize+1));
%    LW = (L-1)*MatrixSize(1);    
%    CentInd = repmat(W',1,length(L))+repmat(LW,length(W),1);
%    CentInd = reshape(CentInd,size(CentInd,1)*size(CentInd,2),1);
end