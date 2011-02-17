function transferContour(imPathNew, imPathOld, imageNameNew, imageNameOld, seriesNameOld, seriesNameNew, indexOld, indexNew, contourName, Trobust)

  
%imNew = imread(imageNameNew);
%imNew = imNew(:,:,1);
%imOld = imread(imageNameOld);
%imOld = imOld(:,:,1);

sizeImOld = [2672,4008];  
sizeImNew = [6219,9329];  

[transformMatrix, ContourPoints, transformMatrixContour, mag, found] = takeTransformFromXmlFile(strcat(imPathOld, seriesNameOld), indexOld, contourName, sizeImOld);

if not(found)
  return;
end

%disp('*********** contour found *********');

for i=1:length(ContourPoints)

cp = ContourPoints(i).points;
cp(:,3) = 100;
cp = cp * inv(transformMatrixContour(i).matrix) * transformMatrix;

scaleFactor = 9329 / 4008;
Tscale = diag([scaleFactor scaleFactor 1]);

cp = cp * Tscale;

cp = cp * Trobust;



fill = ContourPoints(i).fill;
border = ContourPoints(i).border;
mode = ContourPoints(i).mode;

%*******************************
%at this point I have the contour in image coordinates of the newimage
%*********************************

%newPath = 'C:\Dokumente und Einstellungen\Verena\Eigene Dateien\Nuno\Data\4x4 reconstruction\Ubersicht\';

newPath = 'C:\Dokumente und Einstellungen\verena\Eigene Dateien\Nuno\Data\NunoStitched\stitching\';


writeNewContour(newPath, seriesNameNew, indexNew, contourName, cp, sizeImNew, fill, border, mode);
end