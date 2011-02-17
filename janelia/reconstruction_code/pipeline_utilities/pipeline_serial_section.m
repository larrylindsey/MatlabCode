function pipeline_serial_section(config, main_module_id, sub_module_id)
% pipeline_serial_section(module_ids, config)
% Calls different modules of the pipeline based on main_module_id. To be used
% for Serial Section EM Images.
%
% module_ids: 0-999.99
%   0-99 precomputation
%         10  fold detection
%         20  BEL dump training data (TODO)
%         21  BEL dump testing data SAB.
%         22  BEL copy test results into boundary maps (TODO)
%
%   100-199.99 alignment
%         120   SIFT feature point detection
%
%         130   normalized cross corr. tile pair in plane
%         131   normalized cross corr. tile pair inter plane
%
%         150   Deformable mesh correspondence without fold in plane
%         151   SIFT matching tile pair in plane
%         155   Align tile pair in plane
%         158   Display aligned tile pair in plane (TODO)
%
%         160   SIFT matches tile pair inter plane
%         165   Align tile pair inter plane
%         168   Display aligned tile pair inter plane
%
%         175   Joint alignment of all tiles in plane
%         178   Display jointly aligned tiles in plane
%
%         185   Joint alignment of all tiles in 3D stack
%         188   Display jointly aligned tiles in 3D stack
%
%   200-299.99 mitochondria
%         210   manual annotation of regions
%         220   trainig - collect training samples
%         230   train classifier
%         236   train and test in a leave-one-out protocol.
%         240   apply detector
%         246   evaluate detector
%         247   evaluate in a leave-one-out protocol.
%         299   NJS code (TODO).
%
%   400-499.99 superpixel segmentation
%         410     Train segmentation (TODO)
%         450.x   Superpixel sub-module 'x'
%         490     Choose superpixel parameter for each image
%
%   500-599.99 superpixel-to-segment
%         550.x   Suerpixel-to-segent sub-module 'x'
%         590     Choose segmentation parameter for each image
%
%   650   stitch segmentation within a section
%   660   simply assign identity labelling for locked segments from
%           proofread data. No alignment is considered.
%
%   700-799.99 linkage
%         710   linkage train
%         750   linkage apply
%
%   800-899.99  external interface
%         810   export data to proofreader
%           810.1   compute MATLAB style transforms from global alignment
%                     for each section separately.
%           810.15  compile the MATLAB style transforms for all sections.
%           810.2   overlapping area between all the sections.
%           810.3   generate coarse stitched grayscale maps to decide
%                     proofreader ROI.
%           810.35  generate stitched grayscale maps for proofreader
%                     including the ROI.
%           810.39  copy proofread ROI to data to be proofread directory
%                     for records.
%           810.4   generate stitched superpixel maps and
%                     superpixel_2_seg_maps.
%           810.5   combine together the linkage graphs.
%           810.6   post process superpixel_2_seg_map and linkage graph to
%                     ensure that segment ids are unique across stack.
%         850   save reconstruction information.
%         860   launch proofreading session.
%         870   export patch masks to trakEM.
%         880   import data from proofreader
%           880.1   import superpixel map (coarsely) on individual images.
%           880.2   adjust imported superpixel maps in overlapping areas
%                     between tiles.
%           880.21  set the superpixel suffix choice to proofread maps.
%           880.3   import superpixel-2-segment-maps.
%           880.31  set the segmentation suffix choice to proofread data.
%           880.4   import linkage graphs.
%           880.5   import annotations
%           880.55  adjust annotations
%
%   900-999.99 evaluation
%         910   Segment boundary evaluation using bipartite matching.
%           910.1   generate the precision and recall for each pair of
%                     <image, segment suffix> 
%           910.5   compile together the results of 910 into one
%                     precision-recall curve.
%         920   linkage evaluation
%         930   prepare final reconstruction volume with body labels.
%         940   get q score
%         941   get rand score
%         
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% AG mod for linkage
% v0  04162008  init code
% v1  01192009  pipeline version 3
% v2  03112009  branch into 2D and 3D reconstruction
% v3  04212009  major revision of module ids
%

if(config.is_verbose)
  fprintf('START: pipeline_serial_section\n');
end

switch(main_module_id)
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Preprocess images, e.g., fold detection
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(0:99)
    preprocess_images_pipeline(config, main_module_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Alignment
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(100:199)
    alignment_pipeline(config, main_module_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mitochondria - training
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(200:299)
    mitochondria_pipeline(config, main_module_id, sub_module_id);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Superpixel
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(400:499)
    superpixel_main_pipeline(config, main_module_id, sub_module_id);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Superpixel to segment
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(500:599)
    superpixel_2_segment_main_pipeline(config, main_module_id, sub_module_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Alignment of segmentations and images of multiple tiles into
    % one consistent map for each section. This must be performed before
    % linkage.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(600:699)
    segmentation_stich_pipeline(config, main_module_id, sub_module_id)

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linkage
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(700:799)
    linkage_pipeline(config, main_module_id, sub_module_id);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % External interface - export and import data to proofreader, trakEM,
    % etc.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(800:899)
    external_interface_pipeline(config, main_module_id, sub_module_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation modules
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case num2cell(900:999)
    evaluation_pipeline(config, main_module_id, sub_module_id);

  otherwise
    error('Pipeline module id not understood.');
end;

if(config.is_verbose)
  fprintf('STOP: pipeline_serial_section\n');
end
return
end



