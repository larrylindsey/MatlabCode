function import_proofread_data_into_pipeline(config)
% import_proofread_data_into_pipeline(config)
% Imports proofread data from proofreader, e.g., MATLAB gui, Raveler.
%
% Copies onto the image plane the superpixel map from that shown to the
% user during proofreading (possibly stitched from multiple tiles).
% Label '0' is written on all undefined pixels.
%
% Shiv N. Vitaladevuni
% v0  04142009  init. code
%

fprintf('START:import_proofread_data_into_pipeline\n');

switch(config.proofreader.method)
  case 'matlab_gui'
    import_proofread_data_from_mat_file(config);
  case 'Raveler'
    import_proofread_data_from_raveler_v2(config);
  otherwise
    error('Proofreader method not recognized.');
end
fprintf('STOP:import_proofread_data_into_pipeline\n');
return
end
