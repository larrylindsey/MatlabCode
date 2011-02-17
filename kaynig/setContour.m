function setContour(ergStruct, seriesName, indices, imsize, nameOfContours, nonLinear)

  fileName = strcat(seriesName, '.',num2str(indices(1)));
  tree = xmltree(char(fileName));
       
  uidImage = find(tree, '/Section/Transform/Image');
  uidTransform = parent(tree, uidImage);
  mag = attributes(tree,'get',uidImage,1);
  
  mag =str2num(mag.val);
  
  transforms = struct();
  t = diag([1 1 1]);
    
  mag = 0.01 / mag; 
  s = imsize(1)/(100 * mag);
    
  %transform coordinate systems  
		
  tuTrans = [0 mag 0; -mag 0 0; -s 0 1];
  tuRot = [0 1 0; -1 0 0; -s 0 1]; 
    
%  disp('accumulating transforms');
  actualT = diag([1 1 1]);
  ergInd = 0;
  
  for i=1:length(indices)
    disp(i);
    ergInd = ergInd + 1;
%    if or((ergStruct(ergInd).medianErr > 20),(ergStruct(ergInd).matchNumber < 10))
    if (ergStruct(ergInd).matchNumber < 3)
    	 		actualT = diag([1 1 1]);
    	 		disp('not a good match:');
    	 		disp(strcat(['ergStruct: ',num2str(ergInd),' fileExt: ', num2str(indices(i))]))
    else
      actualT = ergStruct(ergInd).Trobust;
    end
    	 
    t = t * actualT;  
    
    tSave = t;
    size = ergStruct(ergInd).imSize;
    sizeDiff = (size(1)-imsize(1)) * 0.01;
    tTransSizeDiff = diag([1 1 1]);
    tTransSizeDiff(3,1) = -sizeDiff;
    
    t = t*tTransSizeDiff;
    
    tRot = t; tRot(3,1) = 0; tRot(3,2) = 0;
    tTrans = diag([1 1 1]); tTrans(3,1) = t(3,1); tTrans(3,2) = t(3,2);
    
    t = tSave;
    
    tNewCoordRot = tuRot * tRot * inv(tuRot);
    
    tNewCoordTrans = tuTrans * tTrans * inv(tuTrans);
    
    
    xcoef = num2str([tNewCoordTrans(3,1)-tNewCoordRot(3,1) tNewCoordRot(1,1) tNewCoordRot(2,1) 0 0 0]);
    ycoef = num2str([tNewCoordTrans(3,2)-tNewCoordRot(3,2) tNewCoordRot(1,2) tNewCoordRot(2,2) 0 0 0]);
    
    index = indices(i);
    fileName = strcat(seriesName, '.',num2str(index));
    fileSaveName = strcat(seriesName, '_',num2str(index),'.old');
    changedTree = xmltree(char(fileName));

    
    
    for contourIndex=1:length(nameOfContours)
      nameOfContour = nameOfContours(contourIndex);
      
      [transformMatrixImageOld, contourPoints, transformMatrixContour, mag, found] = takeTransformFromXmlFile(seriesName, indices(i), nameOfContour, imsize);
      if not(found)
	continue
      end
      for l=1:length(contourPoints)
	cp = contourPoints(l).points;
	tmc = transformMatrixContour(l).matrix;
	cp(:,3) = 100;
	cp = cp*inv(tmc)*transformMatrixImageOld;
	if nonLinear
	  for transInd =1:ergInd
	    [cp, normMean, normVar] = kernelExpandMatrixDimension(cp(:,1:2), 5, ergStruct(transInd).normMean, ergStruct(transInd).normVar, 1);
	    cp = cp * ergStruct(transInd).beta;
	  end
	end
	
	
	cp = cp(:,1:2);
	p1 = (imsize(1)-cp(:,1))*mag;
	p2 = cp(:,2)*mag;
	cp = [p2 p1];
	contour = cp(:,1:2);
	
	
	uidContour = find(changedTree, '/Section/Transform/Contour');
	
	count = 0;
	for j = 1:length(uidContour)
	  name = attributes(changedTree,'get',uidContour(j),1);
	  name = name.val;
	  if strcmp(name,nameOfContour)
	    count = count+1;
	    if count == l
	      points = [];
	      for k=1:length(contour)
		points = strcat([points,num2str(contour(k,1)),' ',num2str(contour(k,2)),',']);
	      end
	      changedTree = attributes(changedTree,'set',uidContour(j),8,'points',points);
	      uidTransform = parent(changedTree, uidContour(j));
	      changedTree = attributes(changedTree,'set',uidTransform,2,'xcoef',xcoef);
	      changedTree = attributes(changedTree,'set',uidTransform,3,'ycoef',ycoef);
	      changedTree = attributes(changedTree,'set',uidTransform,1,'dim','6');
	    end %if count
	  end %if strcmp
	end %for j
      end %for l
    end
    
    uidImage = find(changedTree, '/Section/Transform/Image');
    uidTransform = parent(changedTree, uidImage);
    changedTree = attributes(changedTree,'set',uidTransform,1,'dim','6');
    changedTree = attributes(changedTree,'set',uidTransform,2,'xcoef',xcoef);
    changedTree = attributes(changedTree,'set',uidTransform,3,'ycoef',ycoef);
    
    
    %       save(tree, fileSaveName);
    save(changedTree, fileName);
  end
  