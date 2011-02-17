function deformable_mesh_dir = get_deformable_mesh_dir(config)
% sift_dir = get_deformable_mesh_dir(config)
%
% Get the sift directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  12112008  init code
%

dmesh_config = config.deformable_mesh;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(dmesh_config.dir, 'dir')~=7)
  mkdir2(dmesh_config.dir);
end;
cd(prev_dir);

deformable_mesh_dir = [get_reconstruction_dir(config), dmesh_config.dir];

return
end
