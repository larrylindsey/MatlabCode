%======================================================================
% uidTransform = adjustXMLFilesTraces(seriesName, startIndex, endIndex)
%======================================================================

%takes the reconstruct files according to seriesName and the given
%indexes and changes all transforms found to identitiy.
%uidTransform returns the indexes of all transforms altered.

function uidTransform = adjustXMLFilesTraces(seriesName, startIndex, endIndex)
    for i=startIndex:endIndex
       fileName = strcat(seriesName, '.',num2str(i));
       tree = xmltree(char(fileName));
       
       uidTransform = find(tree, '/Section/Transform');
       changedTree = tree;
       xcoef = num2str([0 1 0 0 0 0]);
       ycoef = num2str([0 0 1 0 0 0]);
       for j=uidTransform
	 changedTree = attributes(changedTree,'set',j,1,'dim','6');
	 changedTree = attributes(changedTree,'set',j,2,'xcoef',xcoef);
	 changedTree = attributes(changedTree,'set',j,3,'ycoef',ycoef);
       end
       
       save(changedTree, fileName);
    end
end
