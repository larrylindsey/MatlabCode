%============================================================
% adjustXMLFiles(ergStruct, seriesName, startIndex, imsize) 
%============================================================
%
% This function takes reconstruct files of the given seriesName and 
% replaces the affine transformation in the files by the
% transformations given in ergStruct.

function adjustXMLFiles(ergStruct, seriesName, startIndex, imsize)

%first calculate accumulated Transforms
  transforms = struct();
  t = diag([1 1 1]);
  
  mag = 0.01 / 0.00254; 
  s = imsize(1)/(100 * mag);
  
  %transform coordinate systems  
  tuTrans = [0 mag 0; -mag 0 0; -s 0 1];
  tuRot = [0 1 0; -1 0 0; -s 0 1]; 
  
  disp('accumulating transforms');
  actualT = diag([1 1 1]);
  
  for i=1:length(ergStruct)
%    if or((ergStruct(i).medianErr > 10),(ergStruct(i).matchNumber <10))
    if false
      actualT = diag([1 1 1]);
      disp('not a good match:');
      disp(i)
    else
      
      scaleDown = 0.25;
      scaleUp = 4.0;

      scaleUp = 1;
      scaleDown = 1;
      actualT = diag([scaleDown scaleDown 1]) *  ergStruct(i).Trobust * diag([scaleUp scaleUp 1]); 
%      actualT = ergStruct(i).Trobust;
    end
    
    t = t * actualT;       
    
    tRot = t; tRot(3,1) = 0; tRot(3,2) = 0;
    tTrans = diag([1 1 1]); tTrans(3,1) = t(3,1); tTrans(3,2) = t(3,2);
       
    tNewCoordRot = tuRot * tRot * inv(tuRot);
       
    tNewCoordTrans = tuTrans * tTrans * inv(tuTrans);
              

    xcoef = num2str([tNewCoordTrans(3,1)-tNewCoordRot(3,1) tNewCoordRot(1,1) tNewCoordRot(2,1) 0 0 0]);
    ycoef = num2str([tNewCoordTrans(3,2)-tNewCoordRot(3,2) tNewCoordRot(1,2) tNewCoordRot(2,2) 0 0 0]);

       
    transforms(i).xcoef = xcoef;
    transforms(i).ycoef = ycoef;
  end
    
  %now read in xml files and adjust parameters
  for i=startIndex:startIndex+length(ergStruct)-1
    fileName = strcat(seriesName, '.',num2str(i));
    fileSaveName = strcat(seriesName, '_',num2str(i),'.old');
    tree = xmltree(char(fileName));
    
% $$$     uidImage = find(tree, '/Section/Transform/Image');
% $$$     uidTransform = parent(tree, uidImage);
% $$$     changedTree = attributes(tree,'set',uidTransform,1,'dim','6');
% $$$     changedTree = attributes(changedTree,'set',uidTransform,2,'xcoef',transforms(i-startIndex+1).xcoef);
% $$$     changedTree = attributes(changedTree,'set',uidTransform,3,'ycoef',transforms(i-startIndex+1).ycoef);

    uidTransform = find(tree, '/Section/Transform');
    
    changedTree = tree;
    
    for uid = 1:length(uidTransform)
      changedTree = attributes(changedTree,'set',uidTransform(uid),1,'dim','6');
      changedTree = attributes(changedTree,'set',uidTransform(uid),2,'xcoef',transforms(i-startIndex+1).xcoef);
      changedTree = attributes(changedTree,'set',uidTransform(uid),3,'ycoef',transforms(i-startIndex+1).ycoef);
    end
    
    %    save(tree, fileSaveName);
    save(changedTree, fileName);
  end
end
