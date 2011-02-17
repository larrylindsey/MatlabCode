function save_reconstruction_param_to_xml(config)
% save_reconstruction_param_to_xml(config)
% save the reconstruction information to xml for records
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

fprintf('Saving reconstruction info. to xml .. ');

xml_save([get_to_be_proofread_dir(config), config.reconstruction.name, '.config.xml'], config);

fprintf('done.\n');

