function [indices, contourNames] = fixXMLFilesAndGetContourNames(seriesName)
  
  dirPath = pwd;
  fileNames = dir(strcat(seriesName,'.*'));

  indices = [];
  contourNames = {};
  
  
  for i=1:length(fileNames)
    name = fileNames(i).name;
    extension = regexprep(name, strcat(seriesName,'.'),'');
    if not(strcmp(extension,'ser'))
      indices = [indices; str2num(extension)];
    end
  end
  
  indices = sort(indices);
  
  for i=1:length(indices)
%    disp(indices(i));
    %for i=find(indices==92)
    fileName = strcat(seriesName,'.',num2str(indices(i)));
    tree = xmltree(char(fileName));
    
    changed = true;
    while (changed)
      changed = false;
      uidTransform = find(tree, '/Section/Transform');
      uidImage = find(tree,'/Section/Transform/Image');
      uidTransformOfImage = parent(tree, uidImage);
    
      if length(uidTransform) > 1
	for j=1:length(uidTransform)
	  if not(uidTransform(j)==uidTransformOfImage)
	    childNodes = children(tree,uidTransform(j));
	    if length(childNodes) > 1
	      tree = copy(tree,uidTransform(j),root(tree));
	      newNode = children(tree,root(tree));
	      newNode = newNode(length(newNode));
	    
	      toDelete1 = childNodes(2:length(childNodes));
	      toDelete2 = children(tree,newNode);
	      toDelete2 = toDelete2(1);
	      
	      tree = delete(tree,[toDelete1 toDelete2]);
	      changed = true;
	    end
	  end
	end
      end
    end %while
    
    
    uidContour = find(tree, '/Section/Transform/Contour');      
    for j=1:length(uidContour)
      name = attributes(tree,'get',uidContour(j),1);
      name = name.val;      
      known = false;
      for k=1:length(contourNames)
	if strcmp(contourNames(k), name) == 1
	  known = true;
	end
      end
      if not(known)
	contourNames(length(contourNames)+1) = {name};
      end
    end
    
    save(tree, fileName);

  end