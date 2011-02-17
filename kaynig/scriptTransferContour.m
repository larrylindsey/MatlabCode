%Needs to be started from old image directory

%imPathNew = 'c:\Dokumente und Einstellungen\Verena\Eigene Dateien\Nuno\Data\4x4 reconstruction\Ubersicht\';

imPathNew = 'C:\Dokumente und Einstellungen\verena\Eigene Dateien\Nuno\Data\NunoStitched\stitching\';

imageNamesNew = dir(strcat(imPathNew,'*.tif'));
length(imageNamesNew)

%imPathOld = 'c:\Dokumente und
%Einstellungen\Verena\Desktop\ContourTransfer\';
imPathOld = 'C:\Dokumente und Einstellungen\verena\Eigene Dateien\Nuno\Data\NunoStitched\contourFiles\';

imageNamesOld = dir(strcat(imPathOld,'*.tif'));
length(imageNamesOld)

seriesNameOld = 'DendriteOverviews';
seriesNameNew = 'test';

[indexOld, contourNames] = fixXMLFilesAndGetContourNames(seriesNameOld);

disp('length indexOld');
length(indexOld)

if not((size(indexOld,1) == length(imageNamesOld)) & (size(indexOld,1) == length(imageNamesNew)))
  disp('ERROR ****** BILDANZAHL NICHT GLEICH!!!!!!!! ***** ERROR');
  return;
end

%ergStruct = struct();

for i=1:length(indexOld)
  
  imageNameNew = strcat(imPathNew, imageNamesNew(i).name)
  imageNameOld = strcat(imPathOld, imageNamesOld(i).name);
  
%  [f1,f2,T, Trobust, iterations, currentErr, s1, sample1, s2, sample2, matchNumber, robustError] = registerHighLow(imageNameNew, imageNameOld, false, 1.0, 4);  

%  ergStruct(i).Trobust = Trobust;
%  ergStruct(i).currentErr = currentErr;
%  ergStruct(i).matchNumber = matchNumber;
  
%  if  (currentErr == -1)
%    continue;
%  end

%Trobust = ergStruct(i).T;%robust;
Trobust = diag([1 1 1]);

  for j=1:length(contourNames)
    cn = contourNames(j);
    transferContour(imPathNew, imPathOld, imageNameNew, imageNameOld, seriesNameOld, seriesNameNew, indexOld(i), i, cn, Trobust);
  end
end
