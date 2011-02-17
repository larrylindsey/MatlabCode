% Doomsday scenario: Manually put correspondence points

image_dir = '/groups/chklovskii/chklovskiilab/electron_microscopy_data/lamina.OTO.shinya.leginon.3500x.July2008/';
dmesh_dir = '../../lamina.OTO.shinya.leginon.3500x.July2008/deformable_mesh/';

% image_sub_dir_1 = 'M2/7-section6/';
% image_prefix_1 = 'M2/7-section6/122_08jul17b_27934M2_6x6_00007gr_00241ex35.mrc';
% image_sub_dir_2 = 'M2/7-section6/';
% image_prefix_2 = 'M2/7-section6/116_08jul17b_27934M2_6x6_00007gr_00247ex35.mrc';

% image_sub_dir_1 = 'M5/5-section11/';
% image_prefix_1 =
% 'M5/5-section11/122_08jul17b_27934M5_6x6_00010gr_00173ex35.mrc';
% image_sub_dir_2 = 'M5/5-section11/';
% image_prefix_2 = 'M5/5-section11/116_08jul17b_27934M5_6x6_00010gr_00179ex35.mrc';

% image_sub_dir_1 = 'N5/6-section6/';
% image_prefix_1 = 'N5/6-section6/128_08jul18a_27934N5_6x6_00007gr_00200ex35.mrc';
% image_sub_dir_2 = 'N5/6-section6/';
% image_prefix_2 = 'N5/6-section6/122_08jul18a_27934N5_6x6_00007gr_00206ex35.mrc';

% image_sub_dir_1 = 'M5/8-section7/';
% image_prefix_1 = 'M5/8-section7/116_08jul17b_27934M5_6x6_00006gr_00287ex35.mrc';
% image_sub_dir_2 = 'M5/8-section7/';
% image_prefix_2 = 'M5/8-section7/115_08jul17b_27934M5_6x6_00006gr_00288ex35.mrc'

% image_sub_dir_1 = 'N5/2-section9/';
% image_prefix_1 =
% 'N5/2-section9/128_08jul18a_27934N5_6x6_00010gr_00056ex35.mrc';
% image_sub_dir_2 = 'N5/2-section9/';
% image_prefix_2 =
% 'N5/2-section9/127_08jul18a_27934N5_6x6_00010gr_00057ex35.mrc';

% image_sub_dir_1 = 'N1/section6/';
% image_prefix_1 = 'N1/section6/115_08jul18a_27934N1_6x6_00007gr_00215ex35.mrc';
% image_sub_dir_2 = 'N1/section6/';
% image_prefix_2 = 'N1/section6/121_08jul18a_27934N1_6x6_00007gr_00209ex35.mrc';
% 
% image_sub_dir_1 = 'P1/4-section6/';
% image_prefix_1 = 'P1/4-section6/122_08jul21a_27934P1_6x6_00007gr_00132ex35.mrc';
% image_sub_dir_2 = 'P1/4-section6/';
% image_prefix_2 = 'P1/4-section6/121_08jul21a_27934P1_6x6_00007gr_00133ex35.mrc';

% image_sub_dir_1 = 'Q2/section4/';
% image_prefix_1 = 'Q2/section4/127_08jul22a_27934Q2_6x6_00007gr_00127ex35.mrc';
% image_sub_dir_2 = 'Q2/section4/';
% image_prefix_2 = 'Q2/section4/121_08jul22a_27934Q2_6x6_00007gr_00133ex35.mrc';
% 

% image_sub_dir_1 = 'Q3/section8/';
% image_prefix_1 = 'Q3/section8/122_08jul22a_27934Q3_6x6_00006gr_00277ex35.mrc';
% image_sub_dir_2 = 'Q3/section8/';
% image_prefix_2 =
% 'Q3/section8/121_08jul22a_27934Q3_6x6_00006gr_00278ex35.mrc';

% image_sub_dir_1 = 'R4/10-section3/';
% image_prefix_1 = 'R4/10-section3/116_08jul23a_27934R4_6x6_00003gr_00355ex35.mrc';
% image_sub_dir_2 = 'R4/10-section3/';
% image_prefix_2 = 'R4/10-section3/115_08jul23a_27934R4_6x6_00003gr_00356ex35.mrc';

% image_sub_dir_1 = 'G4/slice2/';
% image_prefix_1 = 'G4/slice2/121_08jul30a_27934G4_6x6_00003gr_00063ex35.mrc';
% image_sub_dir_2 = 'G4/slice2/';
% image_prefix_2 =
% 'G4/slice2/115_08jul30a_27934G4_6x6_00003gr_00069ex35.mrc';

% image_sub_dir_1 = 'G2/section12/';
% image_prefix_1 = 'G2/section12/133_08jul30a_27934G2_6x6_00010gr_00412ex35.mrc';
% image_sub_dir_2 = 'G2/section12/';
% image_prefix_2 = 'G2/section12/132_08jul30a_27934G2_6x6_00010gr_00413ex35.mrc';
% 

% image_sub_dir_1 = 'F5/slice12/';
% image_prefix_1 = 'F5/slice12/117_08jul31a_27934F5_6x6_00003gr_01057ex35.mrc';
% image_sub_dir_2 = 'F5/slice12/';
% image_prefix_2 = 'F5/slice12/111_08jul31a_27934F5_6x6_00003gr_01063ex35.mrc';

% image_sub_dir_1 = 'F1/slice12/';
% image_prefix_1 = 'F1/slice12/111_08aug01a_27934F1_6x6_00003gr_00435ex35.mrc';
% image_sub_dir_2 = 'F1/slice12/';
% image_prefix_2 = 'F1/slice12/110_08aug01a_27934F1_6x6_00003gr_00436ex35.mrc';
% 
% image_sub_dir_1 = 'E5/slice3/';
% image_prefix_1 = 'E5/slice3/128_08aug01a_27934E5_6x6_00007gr_00093ex35.mrc';
% image_sub_dir_2 = 'E5/slice3/';
% image_prefix_2 = 'E5/slice3/122_08aug01a_27934E5_6x6_00007gr_00099ex35.mrc';


% image_sub_dir_1 = 'O5/section9/';
% image_prefix_1 = 'O5/section9/126_08jul21a_27934O5_6x6_00010gr_00308ex35.mrc';
% image_sub_dir_2 = 'O5/section9/';
% image_prefix_2 = 'O5/section9/120_08jul21a_27934O5_6x6_00010gr_00314ex35.mrc';

image_sub_dir_1 = 'L1/1-section1/';
image_prefix_1 = 'L1/1-section1/121_08jul11b_27934L1_6x6_00004gr_01536ex35.mrc';
image_sub_dir_2 = 'K5/section11/';
image_prefix_2 = 'K5/section11/126_08jul10a_27934K5_6x6_00012gr_00382ex35.mrc';

%%
image1 = im2double(imread([image_dir, image_prefix_1, '.tif']));
image2 = im2double(imread([image_dir, image_prefix_2, '.tif']));

image1 = image1/max(image1(:));
image2 = image2/max(image2(:));

cpselect(image1, image2);
% Two variables are created:
% input_points: image1
% base_points: image2

%%
match_points.points_1 = input_points';
match_points.points_2 = base_points';
file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
  image_prefix_1, image_prefix_2, 'dmcp.');
file_name = [file_name_suffix, '.mat'];
save2(file_name, 'match_points');

match_points.points_1 = base_points';
match_points.points_2 = input_points';
file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
  image_prefix_2, image_prefix_1, 'dmcp.');
file_name = [file_name_suffix, '.mat'];
save2(file_name, 'match_points');
