function external_interface_pipeline(config, main_module_id, sub_module_id)
if(config.is_verbose)
  fprintf('START: external_interface_pipeline\n');
end
switch(main_module_id)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate data structures for proofreading GUI
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 810
    proofread_datastructure_pipeline(config, sub_module_id)

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save reconstruction information
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 850
    save_reconstruction_info(config);
    save_reconstruction_param_to_xml(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Launch proofreading session
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 860
    current_dir = pwd2();
    cd(config.proofreading_session.root_dir);
    switch(sub_module_id)
      case {0, 10}
        start_proofread_session();
      case 20
        start_proofread_session2();
    end
    cd(current_dir);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output datastructures to trakEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 870
    switch(sub_module_id)
      case 1
        output_trakEM_patch_masks_transforms(config);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import data after proofreading
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 880
    switch(sub_module_id)
      case 10
        import_proofread_superpixel_map_from_raveler(config);
      case 20
        adjust_imported_proofread_superpixel_map(config);
      case 21
        set_superpixel_suffix_choice_to_proofread(config);
      case 22
        add_proofread_superpixel_locks(config);
      case 30
        import_proofread_superpixel_2_seg_map(config);
      case 31
        set_segment_suffix_choice_to_proofread(config);
      case 40
        import_proofread_linkage_graph(config);
      case 50
        import_proofread_voxel_annotation_from_raveler(config);
      case 55
        adjust_imported_proofread_voxel_annotation_from_raveler(config);
    end
end
if(config.is_verbose)
  fprintf('STOP: external_interface_pipeline\n');
end
return;
end
