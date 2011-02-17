function IndMat = PatchInd(Center,PatchSize,MatRows)
    % Given a vector of centers and the size of the matrix, return
    % a matrix with the indices for a patch of a given half-size, that is
    % actual size is (2*PatchSize)+1
    
    Base = -PatchSize:PatchSize;
    BaseIndMat = (repmat(Base,length(Base),1)*MatRows)+(repmat(Base,length(Base),1)');
    BaseIndMat = repmat(reshape(BaseIndMat,(2*PatchSize+1).^2,1),1,length(Center));
    IndMat = BaseIndMat+repmat(Center',size(BaseIndMat,1),1);
end