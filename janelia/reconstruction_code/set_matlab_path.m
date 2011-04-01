function set_matlab_path(code_dir)

p = {
  'compare_proofreading', ...
  'segmentation_3D', ...
  'alignment/norm_cross_correlation', ...
  'alignment/tile_align_fold_deformable_mesh', ...
  'alignment/segment_boundary_align', ...
  'alignment', ...
  'binary_io', ...
  'evaluation', ...
  'visualization', ...
  'image_io', ...
  'starting_utilities', ...
  'segmentation_2D', ...
  'segmentation_2D/shrink_segment_levelsets', ...
  'segmentation_2D/boundary_detection', ...
  'segmentation_2D/boundary_detection/LDA', ...
  'segmentation_2D/collect_boundary_stats', ...
  'choose_segmentation_param', ...
  'proofreading_gui', ...
  'preprocess_images', ...
  'pipeline_utilities', ...
  'machine_learning', ...
  'mitochondria_detection', ...
  'misc_utilities', ...
  'linkage_3D', ...
  'raveler_interface_utilities', ...
  'trakEM_interface_utilities', ...
  'house_keeping', ...
  'lib', ...
  'lib/matlab2weka', ...
  'lib/levelsets', ...
  'lib/python_email', ...
  'lib/image_processing_utilities', ...
  'lib/image_processing_matlab', ...
  'lib/graphs_matlab', ...
  'lib/hash_functions', ...
  'lib/misc_matlab_mex', ...
  'lib/misc_matlab_util', ...
  'lib/mouse_draw', ...
  'lib/nearest_neighbor', ...
  'lib/genelib', ...
  'lib/GML_AdaBoost_Matlab_Toolbox_0.3/C++/Example', ...
  'lib/GML_AdaBoost_Matlab_Toolbox_0.3/C++', ...
  'lib/GML_AdaBoost_Matlab_Toolbox_0.3', ...
  'lib/xml_toolbox', ...
  'lib/fftw-3.2', ...
  'lib/fftw-3.2/bin', ...
  'lib/fftw-3.2/include', ...
  'lib/fftw-3.2/lib', ...
  'lib/NcutImage_7_AMD64', ...
  'lib/NcutImage_7_AMD64/common_files', ...
  'lib/NcutImage_7_AMD64/specific_NcutImage_files', ...
  'lib/gaborfilter', ...
  'lib/maxflow-v3.0.src', ...
  'lib/MRF2.1', ...
  'lib/pb', ...
  'lib/readvtk', ...
  'lib/segbench/Benchmark', ...
  'lib/sift', ...
  'lib/parse_json', ...
  };

% code_dir = [strrep(pwd(), '\', '/'), '/'];
for i = 1:length(p)
  dir_name = [code_dir, p{i}];
  if(exist(dir_name, 'dir')==7)
    addpath(dir_name);
  else
    warning('SET_PATH:NONEXISTENT_DIR', ...
      ['Adding path: could not find directory ', dir_name]);
  end
end

javaaddpath([code_dir, 'lib/weka-3-6-2/weka.jar'], '-end');
return
end