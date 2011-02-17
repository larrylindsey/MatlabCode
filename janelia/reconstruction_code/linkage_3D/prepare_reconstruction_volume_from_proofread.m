function recon_volume = prepare_reconstruction_volume_from_proofread()
% recon_volume = prepare_reconstruction_volume_from_proofread()
% Prepare the reconstructed volume after proofreading.
%
% Takes as input cat{}, superpixel_2_seg_map{}, proof and generates a
% volume in which voxel is labelled with a body id.
%
% recon_volume    a cell array of body label maps
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04142008  init code
%

global cat superpixel_2_seg_map proof

n_slice = length(cat);

recon_volume = cell(1, n_slice);

for z = 1:n_slice
  fprintf('%d ', z);
  if(mod(z, 20)==0)
    fprintf('\n');
  end;
  
  % apply the labelling transformations
  recon_volume{z} = double(proof.pmap( superpixel_2_seg_map{z}(cat{z}+1) +1));
  
  % fill in thin 0 label regions within bodies cause due to superpixel
  % boundaries
  recon_volume{z} = imdilate(recon_volume{z}, strel('square', 3));
  
   
  % draw boundaries between superpixels
  boundaries = draw_segment_boundaries_c(recon_volume{z});
  recon_volume{z}(boundaries==1) = 0;
  
%   figure(1); imshow(recon_volume{z}, []);
%   pause;
end



return
end