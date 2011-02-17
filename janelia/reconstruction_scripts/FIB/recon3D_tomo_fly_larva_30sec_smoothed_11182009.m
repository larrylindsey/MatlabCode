function recon3D_tomo_fly_larva_30sec_smoothed_11182009(module_id, case_ids)
% recon3D_fly_larva_FIB_5x5x5nm
%
% module_ids 1000-1999.99
%   1000    precomputation *
%   1001    generate image sub-stack TIFFs after filtering *
%   1002    mitochondria processing *
%   1100    initial 3D watershed *
%   1101    watershed to superpixels *
%   1102    superpixel to segments *
%   1200    generate label TIFF stacks for viewing *
%   1210    output data for proofreading *
%
% * To be implemented
% case_ids    sequence of sections that should be processed (optional)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
% v1  08152008  lamina J recon.
% v2  08252008  test SIFT-based alignment
% v3  10142008  modified for medulla
% v4  12112008  modified to test tile pair alignment under folds using
%                 deformable mesh
% v5  12152008  Example of reconstruction script to be used as template
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

%%% Stack root directory redefinition:
% It is recommended that data be stored in
% chklovskiilab/electron_microscopy_data/. Atypically, the medulla Leginon
% stack is stored in a different directory. This was done because of the
% very large size of the data and the duration of the project.
config.stack.dir = ...
  '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/data_em_images/';
%   '/media/FreeAgent_Drive/em_reconstruction/data_em_images/';
%%% Reconstruction root directory redefinition:
% It is recommended that all intermediate files generated during
% reconstruction be stored in the user's home directory
% em_reconstruction/reconstructions/. However, the directory can redefined
% to a common share if multiple users wish to read and write on the same
% files and to avoid duplicate storage for large projects. For instance,
% uncomment the following line to save intermediate reconstruction files in
% the chklovskii/medulla/ directory:
% config.reconstruction.root_dir = ...
%   '/media/FreeAgent_Drive/em_reconstruction/reconstructions/';
% Directory to store data to be proofread
% config.to_be_proofread.root_dir = ...
%   '/media/FreeAgent_Drive/em_reconstruction/data_to_be_proofread/';
config.DEBUG = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. Configuration parameters - to be set by user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.0 Stack parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) stack name
config.stack.name = 'tomo_fly_larva_30sec_smoothed_11182009';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 's.%03d';  
config.stack.image_suffix = '.tif';
% (c) sections to be analyzed - two options:
% Option 1: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
if(nargin>1)
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 0:84;
end
% % Option 2: List of sub-stacks, specify start and end section numbers.
% config.stack.sub_stack(1).start = 46;
% config.stack.sub_stack(1).end = 75;
% config.stack.sub_stack(2).start = 75;
% config.stack.sub_stack(2).end = 104;
% config.stack.case_ids = ...
%   config.stack.sub_stack(1).start:config.stack.sub_stack(end).end;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi.xmin = 25;
config.stack.roi.xmax = 350;
config.stack.roi.ymin = 25;
config.stack.roi.ymax = 350;
% config.stack.roi = get_roi(config);
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;
% () Set the scale at which segmentation should be performed. This is kept
% constant for all superpixel and superpixel-to-seg stages. When stitching
% the segmentations, the label maps can be further downsampled if needed,
% but not artificially upsampled. During linkage and proofreading, the
% segmentations are upsampled to original resolution.
% config.stack.segmentation_scale = 2;
% config.stack.segmentation_scaling_method = 'bilinear';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.1 Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% name for the reconstruction
config.reconstruction.name = 'test';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1000. Generate image sub-stacks for processing. Useful if full stack is too
% large.
%%%%%%%%%%%%%%%%%%%%%%%%%%
% name for the reconstruction
config.image_stack.filter_version = 'heq_mf3';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1100. Initial 3D watershed 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method to be used for initial segmentation
config.initial_segmentation_3D.method = 'watershed';
% Threshold for watershed (uint8). All saddle-points below this are ignored.
config.initial_segmentation_3D.f_threshold = 70;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1110. Superpixel segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method to be used for superpixel segmentation
config.superpixel_3D.method = 'ladder';
% Method used for generating initial segmentation - see Sec. II.1100
config.superpixel_3D.initial_seg.method = 'watershed';
% Segmentation parameters for initial segmentation - see Sec. II.1100
config.superpixel_3D.initial_seg.suffix = '.ws.T70';
% Threshold on the boundary value when applying ladder (uint8)
config.superpixel_3D.f_threshold = 75;
% Minimum area threshold - segments having area below this value are forced
% to merge.
config.superpixel_3D.area_threshold = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1700. Combine segments from sub-stacks
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Features used for computing adjacency, e.g., area of ovarlap in
% pixels and that normalized by the area in the section.
config.combine_segment_3D.adjacency_method = 'area_overlap_v0';
% many_to_many finds connected components in adjacency graph.
config.combine_segment_3D.link_method = 'many_to_many';
% Parameters for linking overlapping segments
% link = rule_1 | rule_2
% rule_1 = overlap > ta & overlap_norm1 > tb
% rule_2 = overlap_norm2 > tc
% overlap_norm1 = overlap /max(areas)
% overlap_norm2 = overlap /min(areas)
config.combine_segment_3D.area_overlap_threshold = 1000; % ta
config.combine_segment_3D.area_overlap_norm1_threshold = 0.25; % tb
config.combine_segment_3D.area_overlap_norm2_threshold = 0.5; % tc
% Segmentation method to be combined.
config.combine_segment_3D.seg_method = 'watershed';
% Segmentation suffix to be combined.
config.combine_segment_3D.seg_suffix = '.ws.T70';
% Segmentation mapping method to be combined.
config.combine_segment_3D.seg_mapping_method = 'ladder';
% Segmentation mapping suffix to be combined.
config.combine_segment_3D.seg_mapping_suffix = '.ld.T75_L1000';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1900. Display 3D segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whether to print messages
config.display_segmentation_3D.is_verbose = true;
% Segmentation method to be displayed.
config.display_segmentation_3D.seg_method = 'watershed';
% Segmentation suffix to be displayed.
config.display_segmentation_3D.seg_suffix = '.ws.T70';
% % Segmentation mapping method to be displayed.
% config.display_segmentation_3D.seg_mapping_method = 'ladder';
% % Segmentation mapping suffix to be displayed.
% config.display_segmentation_3D.seg_mapping_suffix = '.ld.T75_L1000';
% Segmentation mapping method to be displayed.
config.display_segmentation_3D.seg_mapping_method = 'combined_seg_3D';
% Segmentation mapping suffix to be displayed.
config.display_segmentation_3D.seg_mapping_suffix = '.ld.T75_L1000.aov0_m2m_ta1000_tb0.25_tc0.5';
% Output TIF stack file name.
config.display_segmentation_3D.save_file_name = ['~/temp/test.', ...
  num2str(config.stack.case_ids(1)), '.', ...
  num2str(config.stack.case_ids(end)), '.tif'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end
