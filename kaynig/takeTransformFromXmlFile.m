function [transformMatrix, contourPoints, transformMatrixContour, mag, found] = takeTransformFromXmlFile(seriesName, index, nameOfContour, imsize)

%read in xml files and adjust parameters  
  found = false;
  
  fileName = strcat(seriesName, '.',num2str(index));
  tree = xmltree(char(fileName));
       
  uidImage = find(tree, '/Section/Transform/Image');
  uidTransform = parent(tree, uidImage);
  xcoef = attributes(tree,'get',uidTransform,2);
  ycoef = attributes(tree,'get',uidTransform,3);
  mag = attributes(tree,'get',uidImage,1);
  
  mag =str2num(mag.val);
  
  xcoefContour = [];
  ycoefContour = [];
  
  xcoef = str2num(xcoef.val);
  ycoef = str2num(ycoef.val);
  
  m = 0.01 / mag;
  s = imsize(1)/(100*m);
  
  transformMatrix = diag([1 1 1]);
  transformMatrix(1,1) = ycoef(3);
  transformMatrix(2,2) = xcoef(2);
  transformMatrix(1,2) = -xcoef(3);
  transformMatrix(2,1) = -ycoef(2);

  transformMatrix(3,2) = (xcoef(1)-s*transformMatrix(1,2))*m;
  transformMatrix(3,1) = -(ycoef(1)+s*transformMatrix(1,1)-s)*m;
  
  uidContour = find(tree, '/Section/Transform/Contour');
  transformMatrixContour = struct();
  contourPoints = struct();
  nrContours = 1;
  for i = uidContour
    name = attributes(tree,'get',i,1);
    name = name.val;
    compareLength = (min(length(nameOfContour),length(name)));
    if strcmp(name(1:compareLength),nameOfContour)
      found = true;
      
      mode = attributes(tree,'get',i,7);
      contourPoints(nrContours).mode = mode.val;
      
      fill = attributes(tree,'get',i,6);
      contourPoints(nrContours).fill = fill.val;
      
      border = attributes(tree,'get',i,5);
      contourPoints(nrContours).border = border.val;
      
      points = attributes(tree,'get',i,8);
      
      cp = str2num(points.val);
      uidTransform = parent(tree, i);
      xcoefContour = attributes(tree,'get',uidTransform,2);
      ycoefContour = attributes(tree,'get',uidTransform,3);
    
      xcoefContour = str2num(xcoefContour.val);
      ycoefContour = str2num(ycoefContour.val);
      
      tmc = diag([1 1 1]);
      tmc(1,1) = ycoefContour(3);
      tmc(2,2) = xcoefContour(2);
      tmc(1,2) = -xcoefContour(3);
      tmc(2,1) = -ycoefContour(2);
      
      tmc(3,2) = (xcoefContour(1)-s*tmc(1,2))*m;
      tmc(3,1) = -(ycoefContour(1)+s*tmc(1,1)-s)*m;
      
      transformMatrixContour(nrContours).matrix = tmc;
  
      p1 = cp(1:2:length(cp))/mag;
      p2 = imsize(1) - (cp(2:2:length(cp))/mag);
      cp = [p2;p1]';
      contourPoints(nrContours).points = cp;
      nrContours = nrContours +1;
    end
  end