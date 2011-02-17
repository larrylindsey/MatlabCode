function writeNewContour(path, seriesName, index, nameOfContour, contourPoints, imsize, fill, border, mode)
  
  fileName = strcat(path, seriesName, '.',num2str(index));
  tree = xmltree(char(fileName));
  
  uidImage = find(tree,'/Section/Transform/Image');
  uidTransform = parent(tree, uidImage);

  mag = attributes(tree,'get',uidImage,1);
  mag =str2num(mag.val);

  
  tree = copy(tree,uidTransform,root(tree));

  uidImage = find(tree,'/Section/Transform/Image');
  uidImage = uidImage(length(uidImage));
  uidTransform = parent(tree, uidImage);
  
  tree = delete(tree, uidImage);
  
\..\..\..\Desktop\matlabScriptsForReconstruct\  uidContour = children(tree,uidTransform);
  
  tree = attributes(tree,'set',uidContour,1,'name',nameOfContour);
  tree = attributes(tree,'set',uidContour,2,'hidden','false');
  tree = attributes(tree,'set',uidContour,3,'close','true');
  tree = attributes(tree,'set',uidContour,4,'simplified','true');
  tree = attributes(tree,'set',uidContour,5,'border',border);
  tree = attributes(tree,'set',uidContour,6,'fill',fill);
  tree = attributes(tree,'set',uidContour,7,'mode',mode);

  points = contourPoints;
  
  points(:,1) = imsize(1) - points(:,1) ;
  points = points .* mag;

  
  stringPoints = [];
  for k=1:length(points)
    stringPoints = strcat([stringPoints,num2str(points(k,2)),' ',num2str(points(k,1)),',']);
  end

  
  tree = attributes(tree,'set',uidContour,8,'points',stringPoints);
  save(tree, fileName);