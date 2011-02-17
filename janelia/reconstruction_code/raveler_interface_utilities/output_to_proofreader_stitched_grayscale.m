function output_to_proofreader_stitched_grayscale(al, config)
% output_to_proofreader_stitched_grayscale(al, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

switch(config.proofreader.method)
  case 'matlab_gui'
    save2([get_to_be_proofread_dir(config), 'al.mat'], 'al', '-v7.3');
  case 'Raveler'
    output_to_raveler_stitched_grayscale(al, config)
  otherwise
    error('Proofreader method not recognized');
end

return
end
