function filter_image = boundary_normalized_cuts(image, version_name, config, image_prefix)
% filter_image = boundary_normalized_cuts(image, version_name, config)
% Detect boundaries using normalized cuts
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% Input:
%   image         grayscale image.
%   version_name  specifies the filtering to be performed.
% Outputs:
%   filter_image
%

% perform normalized cuts if not already done
ncut_dir = get_segmentation_ncut_dir(config);
brackets = strfind(version_name, '_');
if(length(brackets)<=2)
  brackets(3) = length(version_name)+1;
end
ncut_version_name = version_name(1:brackets(3)-1);
filter_version_name = version_name(brackets(3)+1:end);

ncut_file_name = ...
  [ncut_dir, image_prefix, '.', ncut_version_name, '.mat'];

if(exist(ncut_file_name, 'file')==2)
  load2(ncut_file_name,  'boundary');
else
  fprintf('\nPerforming normalized cuts - takes a while ...\n');
  scale = config.normalized_cuts.max_image_side/max(size(image));
  if(scale<1)
    width_scaled = round(size(image,2) * scale);
    height_scaled = round(size(image,1) * scale);
    image_s = 255*imresize(image, [height_scaled, width_scaled], 'bilinear');
  else
    image_s = 255*image;
  end

  dataW.sampleRadius=3;
  dataW.sample_rate=1;
  dataW.edgeVariance = 0.1;
  %[number of filter orientations, number of scales, filter size, elongation]
  dataEdgemap.parametres=[4, 3, 15, 5];
  dataEdgemap.threshold=0.002;
  weight_matrix = ICgraph(image_s, dataW, dataEdgemap);

  [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ...
    ncutW(weight_matrix, config.normalized_cuts.n_segments);

  [nr, nc] = size(image_s);
  grad_eigen_vectors = zeros(nr, nc);
  for i=1:config.normalized_cuts.n_segments
    grad_y = NcutEigenvalues(i)*abs(diff(reshape(NcutEigenvectors(:,i) , nr,nc), 1, 1));
    grad_x = NcutEigenvalues(i)*abs(diff(reshape(NcutEigenvectors(:,i) , nr,nc), 1, 2));
    grad_y = [zeros(1, nc); grad_y];
    grad_x = [zeros(nr, 1), grad_x];
    grad = sqrt(grad_x.*grad_x + grad_y.*grad_y);
    grad_eigen_vectors = max(grad_eigen_vectors, grad);
  end

  grad_eigen_vectors = grad_eigen_vectors/max(grad_eigen_vectors(:));
  stack_config = config.stack;
  if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
    image = imresize(image, 1/config.stack.segmentation_scale, ...
      config.stack.segmentation_scaling_method);
  end
  boundary = imresize(grad_eigen_vectors, size(image), 'bilinear'); %#ok<NASGU>
  fprintf('done\n');
  config_normalized_cuts = config.normalized_cuts; %#ok<NASGU>
  save2(ncut_file_name,  'boundary', 'config_normalized_cuts');
end;
filter_image = 1 - boundary;

if(isempty(filter_version_name))
  return;
end
switch(filter_version_name)
  case 'e2'
    filter_image = imerode(filter_image, strel('disk', 2));
  case 'e4'
    filter_image = imerode(filter_image, strel('disk', 4));
  otherwise
    error('Normalized cuts filter option not understood!!!');
end;

return
end

