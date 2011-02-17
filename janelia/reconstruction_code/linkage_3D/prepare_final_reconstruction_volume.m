function prepare_final_reconstruction_volume(config)
% prepare_final_reconstruction_volume(config)
% Prepare the reconstructed volume after proofreading.
%
% Takes as input cat{}, superpixel_2_seg_map{}, proof and generates a
% volume in which voxel is labelled with a body id.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04142008  init code
%

fprintf('Constructing the reconstructed volume ..\n');

seg_config = config.superpixel_2_seg;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;
volume_config = config.final_volume;

load2([get_to_be_proofread_dir(config), 'cat', seg_config.superpixel_suffix, '.mat'], 'cat');
load2([get_to_be_proofread_dir(config), 'superpixel_2_seg_map', apply_config.segmentation_suffix, '.mat'], ...
  'superpixel_2_seg_map');
load2([get_to_be_proofread_dir(config), 'links_3D',  '.',  model_config.type, '.', feature_config.type, ...
  apply_config.model_suffix, apply_config.segmentation_suffix, '.mat'], 'links_3D');

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(volume_config.dir, 'dir')~=7)
  mkdir2(volume_config.dir);
end;
cd(prev_dir);

for t = 1:length(volume_config.link_thresholds)
  fprintf('link threshold = %f\n', volume_config.link_thresholds(t));
  
  pmap = get_pmap_from_linkage_graph(links_3D, volume_config.link_thresholds(t));

  n_slice = length(cat);

  seg = zeros(size(cat{1},1), size(cat{1},2), length(cat));

  for z = 1:n_slice
    fprintf('%d ', z);
    if(mod(z, 20)==0)
      fprintf('\n');
    end;

    % apply the labelling transformations
    recon = double(pmap( superpixel_2_seg_map{z}(cat{z}+1) +1));

%     % fill in thin 0 label regions within bodies cause due to superpixel
%     % boundaries
%     recon = imdilate(recon, strel('square', 3));
% 
% 
%     % draw boundaries between superpixels
%     boundaries = draw_segment_boundaries_c(recon);
%     recon(boundaries==1) = 0;
% 
    seg(:,:,z) = recon;

    %   figure(1); imshow(recon_volume{z}, []);
    %   pause;
  end

  save2([get_reconstruction_dir(config), volume_config.dir, '/', 'seg', '.',  model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, '.lk', num2str(volume_config.link_thresholds(t)), ...
    apply_config.segmentation_suffix, '.mat'], 'seg');

  fprintf('\n');
end;

fprintf('done.\n');


return
end