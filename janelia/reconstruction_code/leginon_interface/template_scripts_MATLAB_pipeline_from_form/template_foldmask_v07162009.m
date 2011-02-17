% Template for calling fold-mask detection module. The user provides
% configuration parameters, e.g., through a web-form, and requests for
% running the fold-mask module. A script based on the following template is
% generated such that it can be run in MATLAB from within
% em_reconstruction/reconstructions/scripts/ directory.
%
% Format of the template:
% - The comments provided in the template are to be replicated in the
% automatically generated scripts.
% - ReconstructionConfigData(<field_name>) is assumed to denote the
% attributeValue for reconstructionName = <reconstruction name> and
% attributeName = <field_name> in table ReconstructionConfigData in the
% reconstruction database. 
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  07162009    first draft
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START: config script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

global code_dir config_global
config_global = [];

%%%%%%%%%%%%%%%%%
% START: Constants not to be changed by user.
%%%%%%%%%%%%%%%%%

%%%%%%%%%
% File name storage - not to be altered.
%%%%%%%%%
% Whether to use hashing when generating file names during storage and
% retrieval. This reduces the filename lengths.
config_global.storage.file_name_is_hashed__ = true;
% Important: Don't change the following values
config_global.storage.hash.string_index_cycle = 4;
config_global.storage.hash.multiplication_factor = (sqrt(5)-1)/2;
config_global.storage.hash.hash_range = 2^20;
% Check to make hash function is compatibile
if(~check_hash_compatibility())
  error('Hash function is not compatible. Saved file names will not be determined.');
end
% During storage, if the file names are to be hashed then part of the file
% name string is hashed. The following define the prefix and suffix of the
% file name that are kept intact for legibility.
% Important: Don't change the following values
config_global.storage.hash.leave_as_is_prefix_length = 50;
config_global.storage.hash.min_length_to_hash = 20;
config_global.storage.hash.leave_as_is_suffix_length = 4;

%%%%%%%%%
% For generating name prefixs from a 2-tuple of strings. Useful for
% generating file names for objects defined over tuples, e.g., transforms
% between pairs of images.  Not to be altered.
%%%%%%%%%
config_global.name_prefix_from_string_tuple.leave_as_is_prefix_length = 6;

%%%%%%%%%
% For deformable mesh based alignment, the valid mesh region labels
% (transform ids) start from TRANSFORMATION_ID_OFFSET. Values 0 to
% TRANSFORMATION_ID_OFFSET are reserved for special regions such as folds,
% etc. Not to be altered.
%%%%%%%%%
config_global.TRANSFORMATION_ID_OFFSET = 10;

%%%%%%%%%
% Standard directories and file name suffixes and extensions.
%%%%%%%%%
% Directory for registering regions.
config.region.registry_dir = ReconstructionConfigData('config.region.registry_dir');

%%%%%%%%%
% Fold mask detection
%%%%%%%%%
% Directory in the stack's reconstruction directory within which the fold
% masks are to be stored.
config.fold.dir = ReconstructionConfigData('config.fold.dir');
% Whether to save the fold masks as TIF files that can be used by
% non-MATLAB programs
config.fold.save_as_tif = true;

%%%%%%%%%%%%%%%%%
% STOP: Constants not to be changed by user.
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%
% START: Global configuration parameters. Affect the behavior of nearly all
% the modules.
%%%%%%%%%%%%%%%%%

%%%%%%%%%
% Whether to execute routines within MATLAB or generate stand alone batch
% scripts that % can be run directly in Linux. Set job.is_stand_along to
% true for generating Linux stand-alone scripts.
%%%%%%%%%
if(~isfield(config_global, 'job') || ...
    ~isfield(config_global.job, 'is_stand_alone'))
  config_global.job.is_stand_alone = false; % default
end

%%%%%%%%%
% Whether running in debug mode. Behavior:
%   print debug messages
%%%%%%%%%
config_global.DEBUG = true;

%%%%%%%%%
% Standard directories and file name suffixes and extensions.
%%%%%%%%%
% Code directory from which programs are to be used.
code_dir = ReconstructionConfigData('code_dir');
% Parent directory within which the recontruction output files are stored
config.reconstruction.parent_dir = ReconstructionConfigData('config.reconstruction.parent_dir');
% Directory within stack's reconstruction directory in which the 2D
% segmentation results are stored
config.superpixel.dir = ReconstructionConfigData('config.segmentation_2D.dir');
config.segmentation_2D.dir = ReconstructionConfigData('config.segmentation_2D.dir');
config.superpixel_2_seg.dir = ReconstructionConfigData('config.segmentation_2D.dir');
% The images are stored in 
% [config.stack.parent_dir, config.stack.name, image_prefix, ...
%   config.stack.image_suffix]
% Parent directory within which the stack directories are stored
config.stack.parent_dir = ReconstructionConfigData('config.stack.parent_dir');
% The extension for images. This is appended to the tile's name.
config.stack.image_suffix = ReconstructionConfigData('config.stack.image_suffix');

%%%%%%%%%
% Reconstruction parameters
%%%%%%%%%
config.reconstruction.name = 'pgrj';

%%%%%%%%%
% Stack parameters
%%%%%%%%%
% name of the stack
config.stack.name = ReconstructionConfigData('config.stack.name');
% trakEM style XML file providing information about the structure of the
% stack: the planes, tilt-planes and the tiles.
config.stack.stack_structure_XML_filename = ReconstructionConfigData('config.stack.stack_structure_XML_filename');

%%%%%%%%%
% Region parameters
%%%%%%%%%
% name of the region
config.region.name = ReconstructionConfigData('config.region.name');
% trakEM style XML file providing information about the structure of the
% region: the planes, tilt-planes and the tiles.
config.region.region_structure_XML_filename = ReconstructionConfigData('config.region.region_structure_XML_filename');
% sections to be processed from the region
config.region.sections = eval(ReconstructionConfigData('config.region.sections'));

%%%%%%%%%
% Region data - structure of the region in terms of plane, tilt-angle and
% tiles.
%%%%%%%%%
config.region.image_data.planes(1).z = 681;
config.region.image_data.planes(1).thickness = 1.0;
config.region.image_data.planes(1).tiltPlanes(1).tilt_angle = 0.0;
config.region.image_data.planes(1).tiltPlanes(1).tiles(1).name = 'box1-M3/6-section1/143_08aug31a_35257M3_9x9_00003gr_00453ex35.mrc.tif';
config.region.image_data.planes(1).tiltPlanes(1).tiles(1).file_path = '/groups/chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/box1-M3/6-section1/143_08aug31a_35257M3_9x9_00003gr_00453ex35.mrc.tif';
config.region.image_data.planes(1).tiltPlanes(1).tiles(2).name = 'box1-M3/6-section1/142_08aug31a_35257M3_9x9_00003gr_00454ex35.mrc.tif';
config.region.image_data.planes(1).tiltPlanes(1).tiles(2).file_path = '/groups/chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/box1-M3/6-section1/142_08aug31a_35257M3_9x9_00003gr_00454ex35.mrc.tif';

config.region.image_data.planes(2).z = 682;
config.region.image_data.planes(2).thickness = 1.0;
config.region.image_data.planes(2).tiltPlanes(1).tilt_angle = 0.0;
config.region.image_data.planes(2).tiltPlanes(1).tiles(1).name = 'box1-M3/7-section2/143_08aug31a_35257M3_9x9_00003gr_00534ex35.mrc.tif';
config.region.image_data.planes(2).tiltPlanes(1).tiles(1).file_path = '/groups/chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/box1-M3/7-section2/143_08aug31a_35257M3_9x9_00003gr_00534ex35.mrc.tif';

%%%%%%%%%
% Fold mask detection
%%%%%%%%%
% Get the configuration from the database's reconstruction config
temp = eval(['struct(', ReconstructionConfigData('config.fold_config'), ')']);
config.fold = cat_struct(config.fold, temp);

%%%%%%%%%
% Superpixel segmentation
%%%%%%%%%
% Get the configuration from the database's reconstruction config
temp = eval(['struct(', ReconstructionConfigData('config.superpixel_config'), ')']);
config.superpixel = cat_struct(config.superpixel, temp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOP: config script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fold mask detection
% main_module_id = 'fold_mask'; sub_module_id = '';
% Superpixel segmentation
% main_module_id = 'superpixel2D'; sub_module_id = '';


pipeline_serial_section(config, main_module_id, sub_module_id);
