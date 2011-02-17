% demoNcutImage
% 
% demo for NcutImage
% also initialize matlab paths to subfolders
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.

clear;
main;

%% read image, change color image to brightness image, resize to 160x160
% I = imread_ncut('specific_NcutImage_files/jpg_images/3.jpg',160,160);
%I = imread_ncut('my_test_images/R34CA1-B_S12.114.crop_2.jpg',320,306);
I = imread_ncut('~/research/segmentation_test_cases/size_450_450.khapical.manualv7.1.jpg',300, 300);

% original_image = im2double(imread('my_test_images/khapical_manualv7_al_55.bmp'));
% n_color = size(original_image, 3);
% if(n_color>1)
%     original_image = rgb2gray(original_image);
% end;
% cropped_image = original_image(1450:1770, 1980:2390);
% [ny nx] = size(cropped_image);
% len = max([ny nx]);
% if(len>400)
%     scale = 400/len;
%     I = imresize(cropped_image, scale, 'bicubic');
% else
%     I = cropped_image;
% end;

%% display the image
figure(1);clf; imagesc(I);colormap(gray);axis off;
disp('This is the input image to segment, press Enter to continue...');
% pause;

%% compute the edges imageEdges, the similarity matrix W based on
% Intervening Contours, the Ncut eigenvectors and discrete segmentation
nbSegments = 50;
disp('computing Ncut eigenvectors ...');
tic;
[SegLabel,NcutDiscrete,NcutEigenvectors,NcutEigenvalues,W,imageEdges]= NcutImage(I,nbSegments);
disp(['The computation took ' num2str(toc) ' seconds on the ' num2str(size(I,1)) 'x' num2str(size(I,2)) ' image']);


%% display the edges
figure(2);clf; imagesc(imageEdges); axis off
disp('This is the edges computed, press Enter to continue...');
pause;

%% display the segmentation
figure(31);clf
bw = edge(SegLabel,0.01);
%J1=showmask(I,imdilate(bw,ones(2,2))); imagesc(J1);axis off
J1=showmask(I,bw); imagesc(J1);axis off
disp('This is the segmentation, press Enter to continue...');
pause;

%% display Ncut eigenvectors
figure(4);clf;set(gcf,'Position',[100,500,200*(nbSegments+1),200]);
[nr,nc,nb] = size(I);
for i=1:nbSegments
    subplot(1,nbSegments,i);
    imagesc(reshape(NcutEigenvectors(:,i) , nr,nc));axis('image');axis off;
end
disp('This is the Ncut eigenvectors...');
disp('The demo is finished.');

