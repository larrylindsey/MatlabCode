% demoNcutImage
%
% demo for NcutImage
% also initialize matlab paths to subfolders
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.

clear;
main;

%% read image, change color image to brightness image, resize to 160x160
% I = imread_ncut('specific_NcutImage_files/jpg_images/3.jpg',160,160);
I = imread_ncut('~/research/cluster_optimization_QSAP/data_videos/horse.1.tif', 160, 160);
% I = rgb2gray(im2double(imread('C:\users\vitaladevunis\n_cuts_experiments\lamina_ex_2.tif')));
% I = imread_ncut('../my_test_images/R34CA1-B_S12.114.crop_2.jpg',320,306);
% I = imread_ncut('../my_test_images/lamina_ex.tif',300,300);
% I = imread_ncut('../my_test_images/lamina_ex.s.tif',263,327);

% maximum_side = 400;
% watershed_init_threshold = 0.001;
% 
% image_file_name = '~/research/normalized_cuts/my_test_images/lamina_ex.tif';
% 
% original_image = imread(image_file_name);
% scale = maximum_side/max(size(original_image));
% width_scaled = round(size(original_image,2) * scale);
% height_scaled = round(size(original_image,1) * scale);
% 
% I = imread_ncut(image_file_name, height_scaled, width_scaled);

[nr, nc] = size(I);

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

%%
dataW.sampleRadius=3;
dataW.sample_rate=1;
dataW.edgeVariance = 0.1;

%[number of filter orientations, number of scales, filter size, elongation]
dataEdgemap.parametres=[4, 3, 6, 3];
dataEdgemap.threshold=0.002;

[W,imageEdges] = ICgraph(I, dataW, dataEdgemap);
figure(2); imshow(imageEdges, []);

%% compute the edges imageEdges, the similarity matrix W based on
% Intervening Contours, the Ncut eigenvectors and discrete segmentation
nbSegments = 60;
disp('computing Ncut eigenvectors ...');
tic;
[SegLabel,NcutDiscrete,NcutEigenvectors,NcutEigenvalues,W,imageEdges]= NcutImage(I,nbSegments);
disp(['The computation took ' num2str(toc) ' seconds on the ' num2str(size(I,1)) 'x' num2str(size(I,2)) ' image']);

% %% display the edges
% figure(2);clf; imagesc(imageEdges); axis off
% disp('This is the edges computed, press Enter to continue...');
% % pause;
% 
% %% display the segmentation
% figure(3);clf
% bw = edge(SegLabel,0.01);
% %J1=showmask(I,imdilate(bw,ones(2,2))); imagesc(J1);axis off
% J1=showmask(I,bw); imagesc(J1);axis off
% disp('This is the segmentation, press Enter to continue...');
% % pause;
% 
% %% display Ncut eigenvectors
% figure(4);clf;set(gcf,'Position',[100,500,200*(nbSegments+1),200]);
% [nr,nc,nb] = size(I);
% for i=1:nbSegments
%   subplot(1,nbSegments,i);
%   imagesc(reshape(NcutEigenvectors(:,i) , nr,nc));axis('image');axis off;
% end
% disp('This is the Ncut eigenvectors...');
% disp('The demo is finished.');

%%
grad_eigen_vectors = zeros(nr, nc);
for i=1:nbSegments
  grad_y = NcutEigenvalues(i)*abs(diff(reshape(NcutEigenvectors(:,i) , nr,nc), 1, 1));
  grad_x = NcutEigenvalues(i)*abs(diff(reshape(NcutEigenvectors(:,i) , nr,nc), 1, 2));
  grad_y = [zeros(1, nc); grad_y];
  grad_x = [zeros(nr, 1), grad_x];
  grad = sqrt(grad_x.*grad_x + grad_y.*grad_y);
  grad_eigen_vectors = max(grad_eigen_vectors, grad);
end

grad_eigen_vectors = grad_eigen_vectors/max(grad_eigen_vectors(:));
grad_eigen_vectors = imresize(grad_eigen_vectors, size(original_image), 'bilinear');
figure(5); imshow(grad_eigen_vectors, []);

%%
grad_thresh = grad_eigen_vectors.*(grad_eigen_vectors>watershed_init_threshold);
g_ws = watershed(grad_thresh, 4);
figure(6); imshow(g_ws==0);

%%
f_threshold = 0.01;
area_threshold = 150;
l = compute_segmentation_hierarchy_from_watershed_with_min_area_4c(...
  g_ws, grad_eigen_vectors, f_threshold, area_threshold);
figure(7); imshow(original_image, []);
[py, px] = find(l==0);
hold on; plot(px, py, '.'); hold off;
