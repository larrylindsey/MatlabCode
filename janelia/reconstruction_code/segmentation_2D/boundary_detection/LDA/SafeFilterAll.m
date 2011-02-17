function [Val,CentInd] = SafeFilterAll(Matrix,PatchSize,Filter)

fprintf('START: SafeFilterAll\n');
% Apply a filter to all pixels in the image
OMS = size(Matrix);
%   OM = Matrix;
Matrix = padarray(Matrix,[PatchSize PatchSize],'symmetric');
display('Calculating all pixel filter values')
FilterOutputSize = size(Filter,2);
MatrixSize = size(Matrix);
%    W = (PatchSize+1):(MatrixSize(1)-(PatchSize+1));
%    L = (PatchSize+1):(MatrixSize(2)-(PatchSize+1));
% 25.2.2010
W = (PatchSize+1):(MatrixSize(1)-(PatchSize));
L = (PatchSize+1):(MatrixSize(2)-(PatchSize));

LW = (L-1)*MatrixSize(1);
%    WW = (W-1)*MatrixSize(1);
%    CentInd = repmat(W',1,length(W))+repmat(WW,length(W),1);
CentInd = repmat(W',1,length(L))+repmat(LW,length(W),1);
CentInd = reshape(CentInd,size(CentInd,1)*size(CentInd,2),1);
[I,J] = ind2sub(size(Matrix),CentInd);
I = I-PatchSize;J = J-PatchSize;
NewCentInd = sub2ind(OMS,I,J);
Val = zeros(length(CentInd),FilterOutputSize);

PatchArea = (2*PatchSize+1).^2;
MemLim = 1e7;
PatchLim = floor(MemLim/PatchArea);
if length(CentInd) < PatchLim
  PatchSequenceVec = [1 length(CentInd)];
else
  PatchSequenceVec = [1:PatchLim:length(CentInd) length(CentInd)];
  if(PatchSequenceVec(end-1) == PatchSequenceVec(end))
    PatchSequenceVec(end) = [];
  end
end

for ii=1:(length(PatchSequenceVec)-1)
  CurInd = PatchSequenceVec(ii):PatchSequenceVec(ii+1);
  Patch = Matrix(PatchInd(CentInd(CurInd),PatchSize,MatrixSize(1)));
  %       Val(CurInd,:) = FilterHandle(Patch',Filter);
  Val(CurInd,:) = Patch'*Filter;
  
end
% Seems like the padding works
CentInd = NewCentInd;

fprintf('STOP: SafeFilterAll\n');

return;
end
