function copy_global_patchwise_transforms_from_simple_txt_to_mat_files(...
  region_structure_file_name)
fprintf('START: copy_global_patchwise_transforms_from_simple_txt_to_mat_files\n');
config = get_basic_config();

input_file_name = '~/temp/simple.02262010';
input_fold_mask_dir = '/groups/scheffer/home/schefferl/Pictures/temp/stack/';
tform_output_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/region.crop4_global_alignment_0161_1460.161.1460.c352854878051_40/deformable_mesh.02262010/';
fold_mask_output_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/fold_masks.02262010/';

config.stack.dir = '/groups/chklovskii/medulla/';
config.stack.name = 'medulla.HPF.Leginon.3500x.zhiyuan.fall2008';   
config.stack.image_structure = 'crop4_global_alignment_0161_1460.kludge.xml';
config.stack.image_structure_is_kludged = true;

config.region.region_structure = load(region_structure_file_name);


for i = 1:length(config.region.region_structure.planes)
  tilt_plane = config.region.region_structure.planes(i).tilt_planes(1);
  grep_string = tilt_plane.tiles(1).sub_dir(2:end);
  for j = 1:length(tilt_plane.tiles)
    tile = tilt_plane.tiles(j);
    fprintf('tile: %s\n', tile.name);
    
    comm = ['grep ', grep_string, ' ', input_file_name, ' | grep "FOLDMAP" | grep "', tile.name, '" | uniq | ', ...
      'awk ''{print $3;}'''];
    [status, result] = system(comm);
    if(status~=0)
      error('Error while recovering fold mask file name');
    end
    result = result(1:end-1);
    fold_mask_filename = [input_fold_mask_dir, result];
    fold_mask_filename = [fold_mask_filename(1:end-4), 'd.tif'];
    fprintf('fold_mask_filename: %s\n', fold_mask_filename);
    fold_mask = imread(fold_mask_filename);
    
    n_patch = max(fold_mask(:));
    fprintf('n_patch: %d\n', n_patch);
    
    fold_mask_new = zeros(size(fold_mask));
    for patch_id = 1:n_patch
      mask = fold_mask==patch_id;
      mask = imfill(mask, 4, 'holes');
      fold_mask_new = max(fold_mask_new, double(mask)*double(patch_id));
    end
    fold_mask = uint8(fold_mask_new);
    
    check_for_dir([fold_mask_output_dir, tile.sub_dir]);
%     fold_mask_output_file_name = ...
%       [fold_mask_output_dir, tile.name, '.fold_mask.mat'];
%     save2(fold_mask_output_file_name, 'fold_mask');
    fold_mask_output_file_name = ...
      [fold_mask_output_dir, tile.name, '.fold_mask.tif'];
    imwrite(fold_mask, get_storage_file_name(fold_mask_output_file_name));
    
    transform.transform = [];
    for patch_id = 1:n_patch
      comm = ['grep ', grep_string, ' ', input_file_name, ' | grep "TRANSFORM" | grep "', tile.name, '.tif"\''"::', ...
        num2str(patch_id), '" | uniq | ', ... 
        'awk ''{print $3 " " $4 " " $5 " " $6 " " $7 " " $8}'''];
      [status, result] = system(comm);
      if(status~=0)
        error('Error while recovering transform');
      end
      result = result(1:end-1);
      if(isempty(result))
        transform.transform(:,end+1) = zeros(6, 1);
      else
        transform.transform(:,end+1) = str2num(result)'; %#ok<ST2NM>
        transform.transform(:,end) = transform.transform([1 4 2 5 3 6],end);
      end
    end
    fprintf('transform:\n');
    display(transform);
    is_empty_tform = max(abs(transform.transform))==0;
    if(nnz(is_empty_tform)>0)
      copy_patch = find(~is_empty_tform, 1);
      transform.transform(:, is_empty_tform) = ...
        repmat(transform.transform(:, copy_patch), [1, nnz(is_empty_tform)]);
      fprintf('transform (filled in) :\n');
      display(transform.transform);
    end
    
    check_for_dir([tform_output_dir, tile.sub_dir]);
    tform_output_file_name = ...
      [tform_output_dir, tile.name, '.patchwise_affine_transforms.mat'];
    save2(tform_output_file_name, 'transform');
  end
end

fprintf('STOP: copy_global_patchwise_transforms_from_simple_txt_to_mat_files\n');
return
end
