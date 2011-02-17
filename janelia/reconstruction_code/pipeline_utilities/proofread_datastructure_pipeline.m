function proofread_datastructure_pipeline(config, sub_module_id)
if(config.is_verbose)
  fprintf('START: proofread_datastructure_pipeline\n');
end
switch(sub_module_id)
  case 0
    if(isfield(config.stack, 'image_structure') ...
        && ~isempty(config.stack.image_structure))
      generate_matlab_global_tforms(config);
      generate_matlab_global_tforms_compile(config);
      generate_align_roi(config);
      generate_al(config);
      copy_proofread_roi_to_output_dir(config);
      generate_cat_superpixel_2_seg_map(config);
      combine_linkage_graphs_multi_tile(config);
      process_superpixel_2_2Dseg_linkage_graph(config);
    else
      generate_matlab_global_tforms(config);
      generate_matlab_global_tforms_compile(config);
      generate_align_roi(config);
      generate_al(config);
      copy_proofread_roi_to_output_dir(config);
      generate_cat_superpixel_2_seg_map(config);
      combine_linkage_graphs_uni_tile(config);
      process_superpixel_2_2Dseg_linkage_graph(config);
    end
  case 5
    dump_info_for_c_rendering(config);
  case 10
    generate_matlab_global_tforms(config);
  case 15
    generate_matlab_global_tforms_compile(config);
  case 20
    generate_align_roi(config);
    
  case 30
    generate_coarse_al_without_roi(config);
  case 35
    generate_al(config);
  case 37
    generate_annotations(config);
  case 39
    copy_proofread_roi_to_output_dir(config);
    
  case 40
    generate_cat_superpixel_2_seg_map(config);
  case 50
    combine_linkage_graphs(config);
  case 60
    process_superpixel_2_2Dseg_linkage_graph(config);
end
if(config.is_verbose)
  fprintf('STOP: proofread_datastructure_pipeline\n');
end
return
end
