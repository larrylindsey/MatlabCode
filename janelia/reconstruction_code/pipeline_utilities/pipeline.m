function pipeline(module_ids, config)
% pipeline(module_ids, config)
% Calls different modules of the pipeline based on module_id
%
% module_id
%   0-999.99        For reconstruction of Serial Section Electron
%                     Micrographs. Type help pipeline_serial_section.
%   999.99-1999.99  For reconstruction of Block Face and Focus Ion Beam
%                     Micrographs. Type help pipeline_block_face.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% AG mod for linkage
% v0  04162008  init code
% v1  01192009  pipeline version 3
% v2  03112009  branch into 2D and 3D reconstruction
%

config = process_config_structure(config);

prepare_job(config);

module_ids = int32(module_ids*100);

for module_id = module_ids
  main_module_id = int32(floor(double(module_id)/100));
  sub_module_id = int32(mod(module_id,100));

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Pipeline for Serial Section Electron Micrographs
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  if(main_module_id>=0 && main_module_id<999)
    pipeline_serial_section(config, main_module_id, sub_module_id);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Pipeline for FIB or Block Face Electron Micrographs
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  if(main_module_id>=1000 && main_module_id<2000)
    pipeline_block_face(config, main_module_id, sub_module_id);
  end
end;

pipeline_cleanup(config);

beep
return
end
