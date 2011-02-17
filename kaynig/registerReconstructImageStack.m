function ergsStruct = registerReconstructImageStack(dirPath, seriesName, startIndex)
imageNames = dir(strcat(dirPath,'*.tif'));

ergsStruct = struct();

ergsStruct(1).T = diag([1 1 1]);
ergsStruct(1).Trobust = diag([1 1 1]);

frames1 = [];
descriptors1 = [];

for i=2:length(imageNames);
  disp(imageNames(i-1).name);
  imPath1 = strcat(dirPath,imageNames(i-1).name);
  imPath2 = strcat(dirPath,imageNames(i).name);
  
  [T, Trobust, iterations, currentErr, sample1, sample2, matchNumber, robustError, frames1, descriptors1] = registrationSIFTRigid(imPath1, imPath2, false,frames1, descriptors1);
  
    
  ergsStruct(i).T = Trobust;
  ergsStruct(i).Trobust = T;
%  ergsStruct(i).iterations = iterations;
  ergsStruct(i).medianErr = currentErr;
%vor huge data sets saving the sample points can lead to out of memory
%messages
%  ergsStruct(i).sample1 = sample1;
%  ergsStruct(i).sample2 = sample2;
  ergsStruct(i).matchNumber = matchNumber;
  ergsStruct(i).robustError = robustError;

  save rigidRegistration.mat ergsStruct;
end

tmpIm = imread(strcat(dirPath,imageNames(1).name));
imsize = size(tmpIm);

disp('Adjusting XML Files...');

adjustXMLFiles(ergsStruct, seriesName, startIndex, imsize) 



