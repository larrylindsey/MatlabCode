function deformable_mesh_dir = get_norm_x_corr_dir(config)
% sift_dir = get_deformable_mesh_dir(config)
%
% Get the sift directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  12112008  init code
%

norm_x_corr_config = config.norm_x_corr;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(norm_x_corr_config.dir, 'dir')~=7)
  mkdir2(norm_x_corr_config.dir);
end;
cd(prev_dir);

deformable_mesh_dir = [get_reconstruction_dir(config), norm_x_corr_config.dir];

return
end
