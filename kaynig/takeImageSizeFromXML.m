function [size, found] = takeImageSizeFromXML(seriesName, index, nameOfContour)

%read in xml files and adjust parameters  
  found = false;
  
  fileName = strcat(seriesName, '.',num2str(index));
  tree = xmltree(char(fileName));
       
  uidImage = find(tree, '/Section/Transform/Image');
  uidTransform = parent(tree, uidImage);
  
  uidContour = find(tree, '/Section/Transform/Contour');
  nrContours = 1;
  for i = uidContour
    name = attributes(tree,'get',i,1);
    name = name.val;
    if strcmp(name,nameOfContour)
      found = true;
      points = attributes(tree,'get',i,8);
      cp = str2num(points.val);

      p1 = cp(1:2:length(cp));
      p2 = cp(2:2:length(cp));
      cp = [p2;p1]';
      size = max(cp);
    end
  end