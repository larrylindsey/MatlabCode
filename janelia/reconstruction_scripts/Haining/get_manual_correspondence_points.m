% Doomsday scenario: Manually put correspondence points

image_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/data_em_images/Haining.08142009/';
alignment_dir = '../../Haining.08142009/alignmentSIFT/';

image_prefix_1 = '08142009/images/Fig1';
image_prefix_2 = '08142009/images/Fig2';

%%
image1 = im2double(imread([image_dir, image_prefix_1, '.tif']));
image2 = im2double(imread([image_dir, image_prefix_2, '.tif']));

image1 = image1*10;
image2 = image2*10;

cpselect(image1, image2);
% Two variables are created:
% input_points: image1
% base_points: image2

%%
match_points.points_1 = input_points';
match_points.points_2 = base_points';
file_name_suffix = [alignment_dir, 'match_points.', ...
  image_prefix_1, '.', image_prefix_2];
file_name = [file_name_suffix, '.mat'];
save2(file_name, 'match_points');

match_points.points_1 = base_points';
match_points.points_2 = input_points';
file_name_suffix = [alignment_dir, 'match_points.', ...
  image_prefix_2, '.', image_prefix_1];
file_name = [file_name_suffix, '.mat'];
save2(file_name, 'match_points');

%%
p1 = input_points';
p1 = [p1; ones(1, size(p1,2))];

p2 = base_points';
p2 = [p2; ones(1, size(p2,2))];

A = p1*p2'/(p2*p2');

transform = [1 0 0 1 0 0]; %#ok<NASGU>
file_name = [alignment_dir, image_prefix_1, '.affine_transforms.mat'];
save2(file_name, 'transform');

transform = reshape(A(1:2,:), [1 6]);
file_name = [alignment_dir, image_prefix_2, '.affine_transforms.mat'];
save2(file_name, 'transform');

