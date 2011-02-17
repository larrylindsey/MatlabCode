function output_trakEM_patch_masks_transforms(config)
% output_trakEM_patch_masks_transforms(config)
% Output patch masks and their global transforms to trakEM. Generates a
% javascript that can be run within trakEM to manipulate a project.
%
%
% Input:
%   config    config datastructure of the reconstruction
%
% trakEM interface provided by Stephan Saalfeld and Albert Cardona.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  01022009  init. code
%

stack_config = config.stack;
trakEM_config = config.trakEM.output_patch_masks_transforms;

image_dir = get_stack_dir(config);
fold_dir = [get_reconstruction_dir(config), trakEM_config.dir];
switch(config.align.global_align.method)
  case 'SIFT'
    tform_dir = get_sift_dir(config);
  case 'deformable_mesh'
    tform_dir = get_deformable_mesh_dir(config);
end

if( ~isfield(stack_config, 'image_structure') || isempty(stack_config.image_structure) )
  error('xml file from TrakEM not specified. Could not load tile images and initial transformations.');
end;

if(~isfield(trakEM_config, 'is_verbose'))
  trakEM_config.is_verbose = true;
end
if(trakEM_config.is_verbose)
  fprintf('Computing fold masks ..\n');
end

s = xmlread([image_dir, stack_config.image_structure]);
xmlwrite('~tmp.xml', s);
xmlstr = fileread('~tmp.xml');
s = xml_parseany(xmlstr);
s_t2_layer = s.t2_layer_set{1}.t2_layer;

if(exist(trakEM_config.output_script_name, 'file')==2)
  delete(trakEM_config.output_script_name);
end
copyfile(trakEM_config.preamble_file, trakEM_config.output_script_name, 'f');
fout = fopen(trakEM_config.output_script_name, 'at');

for i = 1:length(s_t2_layer)
  case_id=str2double(s_t2_layer{i}.ATTRIBUTE.z);
  if(ismember(case_id,stack_config.case_ids) || isempty(stack_config.case_ids))
    if(trakEM_config.is_verbose)
      fprintf('%d\n', case_id);
    end
    for j = 1:length(s_t2_layer{i}.t2_patch)
      % remove the directory and extension from the file_path and store
      % in image_prefix
      dir_length = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '/');
      dir_length = [zeros(1, 3-length(dir_length)), dir_length];

      % File path in .xml file is
      % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
      image_prefix = ...
        s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end-4);

      file_name_prefix = [fold_dir, image_prefix, '.fold_mask'];
      fold_mask = load_fold_mask(file_name_prefix, config);
      n_patch = max(fold_mask(:));
      
      load2([tform_dir, image_prefix, '.patchwise_affine_transforms.mat'], 'transform');
      
      if(n_patch==1)
        % entire image is one patch - just update the transformation
        fprintf(fout, '\n');
        fprintf(fout, ['\ttile = tiles.get( "', ...
          s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '" );\n']);
        fprintf(fout, '\ttile.setAffineTransform( new AffineTransform( ');
        for param_id = 1:5
          fprintf(fout, '%g, ', transform.transform(param_id));
        end
        fprintf(fout, '%g ) );\n', transform.transform(6));
        fprintf(fout, '\ttile.updateMipmaps();\n');
      else
        % save patch files for trakEM
        for patch_id = 1:n_patch
          patch_mask = 255*(fold_mask == patch_id);
          file_name = [fold_dir, image_prefix, '.patch_mask.', num2str(patch_id), '.tif'];
          imwrite(patch_mask, file_name);
        end
        % multiple patches - insert patch and transform info.
        % first patch
        fprintf(fout, '\n');
        patch_id = 1;
        fprintf(fout, ['\timpMask = opener.openImage( "', ...
          fold_dir, image_prefix, '.patch_mask.', num2str(patch_id), '.tif" );\n']);
        fprintf(fout, ['\ttile = tiles.get( "', ...
          s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '" );\n']);
        fprintf(fout, '\ttile.setAffineTransform( new AffineTransform( ');
        for param_id = 1:5
          fprintf(fout, '%g, ', transform.transform(param_id, patch_id));
        end
        fprintf(fout, '%g ) );\n', transform.transform(6, patch_id));
        fprintf(fout, '\ttile.setAlphaMask( impMask.getProcessor() );\n');
        fprintf(fout, '\ttile.updateMipmaps();\n');
        % rest of the patches
        for patch_id = 2:n_patch
          fprintf(fout, ['\timpMask = opener.openImage( "', ...
            fold_dir, image_prefix, '.patch_mask.', num2str(patch_id), '.tif" );\n']);
          fprintf(fout, '\ttileCopy = tile.clone();\n');
          fprintf(fout, '\ttile.getLayer().add( tileCopy );\n');
          fprintf(fout, '\ttileCopy.setAffineTransform( new AffineTransform( ');
          for param_id = 1:5
            fprintf(fout, '%g, ', transform.transform(param_id, patch_id));
          end
          fprintf(fout, '%g ) );\n', transform.transform(6, patch_id));
          fprintf(fout, '\ttileCopy.setAlphaMask( impMask.getProcessor() );\n');
          fprintf(fout, '\ttileCopy.updateMipmaps();\n');
        end        
      end
    end
  end
end

fin = fopen(trakEM_config.postscript_file, 'rt');
fprintf(fout, '%c', fscanf(fin, '%c', inf));
fclose(fin);

fclose(fout);

if(trakEM_config.is_verbose)
  fprintf('done.\n');
end

return
end
