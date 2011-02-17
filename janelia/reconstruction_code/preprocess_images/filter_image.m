function filtered_image = filter_image(image, version_name, config, image_prefix)
% [filtered_image] = filtered_image = filtered_image(image, version_name, config,
%                   image_prefix)
% Filter the image to reduce noise prior to segmentation
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% Input:
%   image         grayscale image.
%   version_name  specifies the filtering to be performed.
%   config        reconstruction configuration (optional)
%   image_prefix  image's file name prefix for saving intermediate results,
%                   e.g., normalized cuts. (optional).
% Outputs:
%   filtered_image
%

global config_global
if(isfield(config_global, 'DEBUG') && config_global.DEBUG)
  is_verbose = true;
else
  is_verbose = false;
end

if(is_verbose)
  fprintf('\nSTART: filter_image\n');
end

if(~iscell(version_name))
  version_names = {version_name};
else
  version_names = version_name;
end

for v_id = 1:length(version_names)
  version_name = version_names{v_id};
  if(is_verbose)
    fprintf('version_name: %s\n', version_name);
  end

  if(strcmp(version_name(1), 'v')~=1)
    if(nargin==4)
      filtered_image = filter_image2(image, version_name, config, image_prefix);
    elseif(nargin==3)
      filtered_image = filter_image2(image, version_name, config);
    else
      filtered_image = filter_image2(image, version_name);
    end
    image = filtered_image;
    continue;
  end
  
  % For backward compatibility
  switch(version_name)
    case 'v0'
      filtered_image = image;
    case 'v_heq2'
      filtered_image = histeq2(image(image>0.1 & image<0.9), image);
    case 'v_heq'
      filtered_image = histeq(image);
    case 'v_heq2_mf2'
      image = histeq2(image(image>0.1 & image<0.9), image);
      filtered_image = medfilt2(image, [2 2], 'symmetric');
    case 'v_heq2_mf2_e2'
      image = histeq2(image(image>0.1 & image<0.9), image);
      image = medfilt2(image, [2 2], 'symmetric');
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_heq2_mf2_e1'
      image = histeq2(image(image>0.1 & image<0.9), image);
      image = medfilt2(image, [2 2], 'symmetric');
      filtered_image = imerode(image, strel('disk', 1));
    case 'v_heq2_mf10'
      image = histeq2(image(image>0.1 & image<0.9), image);
      filtered_image = medfilt2(image, [10 10], 'symmetric');
    case 'v_aheq'
      filtered_image = adapthisteq(image);
    case 'v_heq_mf1_e2'
      image = histeq(image);
      image = medfilt2(image, [1 1], 'symmetric');
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_heq_mf1'
      image = histeq(image);
      filtered_image = medfilt2(image, [1 1], 'symmetric');
    case 'v_heq_mf2'
      image = histeq(image);
      filtered_image = medfilt2(image, [2 2], 'symmetric');
    case 'v_heq_mf2_e1'
      image = histeq(image);
      image = medfilt2(image, [2 2], 'symmetric');
      filtered_image = imerode(image, strel('disk', 1));
    case 'v_heq_mf2_e2'
      image = histeq(image);
      image = medfilt2(image, [2 2], 'symmetric');
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_heq_mf3'
      image = histeq(image);
      filtered_image = medfilt2(image, [3 3], 'symmetric');
    case 'v_heq_mf4'
      image = histeq(image);
      filtered_image = medfilt2(image, [4 4], 'symmetric');
    case 'v_heq_mf5'
      image = histeq(image);
      filtered_image = medfilt2(image, [5 5], 'symmetric');
    case 'v_heq_mf10'
      image = histeq(image);
      filtered_image = medfilt2(image, [10 10], 'symmetric');
    case 'v_heq_mf3_e1'
      image = histeq(image);
      image = medfilt2(image, [3 3], 'symmetric');
      filtered_image = imerode(image, strel('disk', 1));
    case 'v_heq_mf3_e2'
      image = histeq(image);
      image = medfilt2(image, [3 3], 'symmetric');
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_heq_mf3_e3'
      image = histeq(image);
      image = medfilt2(image, [3 3], 'symmetric');
      filtered_image = imerode(image, strel('disk', 3));
    case 'v_heq_mf4_e2'
      image = histeq(image);
      image = medfilt2(image, [4 4], 'symmetric');
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_heq_mf5_e1'
      image = histeq(image);
      image = medfilt2(image, [5 5], 'symmetric');
      filtered_image = imerode(image, strel('disk', 1));
    case 'v_heq_mf5_e2'
      image = histeq(image);
      image = medfilt2(image, [5 5], 'symmetric');
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_heq_mf5_e3'
      image = histeq(image);
      image = medfilt2(image, [5 5], 'symmetric');
      filtered_image = imerode(image, strel('disk', 3));
    case 'v_heq_mf5_d6'
      image = histeq(image);
      image = medfilt2(image, [5 5], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 6));
    case 'v_heq_mf10_e1'
      image = histeq(image);
      image = medfilt2(image, [10 10], 'symmetric');
      filtered_image = imerode(image, strel('disk', 1));
    case 'v_heq_mf10_e2'
      image = histeq(image);
      image = medfilt2(image, [10 10], 'symmetric');
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_heq_mf10_e3'
      image = histeq(image);
      image = medfilt2(image, [10 10], 'symmetric');
      filtered_image = imerode(image, strel('disk', 3));
    case 'v_heq_e1'
      image = histeq(image);
      filtered_image = imerode(image, strel('disk', 1));
    case 'v_heq_e2'
      image = histeq(image);
      filtered_image = imerode(image, strel('disk', 2));
    case 'v_mf2' % for boundary detectors such as BEL, etc.
      filtered_image = medfilt2(image, [2 2], 'symmetric');
    case 'v_mf3' % for boundary detectors such as BEL, etc.
      filtered_image = medfilt2(image, [3 3], 'symmetric');
    case 'v_mf4' % for boundary detectors such as BEL, etc.
      filtered_image = medfilt2(image, [4 4], 'symmetric');
    case 'v_mf3_d1' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [3 3], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 1));
    case 'v_mf3_d2' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [3 3], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 2));
    case 'v_mf5' % for boundary detectors such as BEL, etc.
      filtered_image = medfilt2(image, [5 5], 'symmetric');
    case 'v_mf5_d1' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [5 5], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 1));
    case 'v_mf5_d2' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [5 5], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 2));
    case 'v_mf5_d3' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [5 5], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 3));
    case 'v_mf10' % for boundary detectors such as BEL, etc.
      filtered_image = medfilt2(image, [10 10], 'symmetric');
    case 'v_mf10_d1' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [10 10], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 1));
    case 'v_mf10_d2' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [10 10], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 2));
    case 'v_mf10_d3' % for boundary detectors such as BEL, etc.
      image = medfilt2(image, [10 10], 'symmetric');
      filtered_image = imdilate(image, strel('disk', 3));
    case 'v_mf15' % for boundary detectors such as BEL, etc.
      filtered_image = medfilt2(image, [15 15], 'symmetric');
    case 'v_heq_o1'
      image = histeq(image);
      filtered_image = imopen(image, strel('disk', 1));
    case 'v_heq_d1'
      image = histeq(image);
      filtered_image = imdilate(image, strel('disk', 1));
    case 'v_heq_d2'
      image = histeq(image);
      filtered_image = imdilate(image, strel('disk', 2));
    case 'v_heq_mf2_d1'
      image = histeq(image);
      filtered_image = imdilate(medfilt2(image, [2 2]), strel('disk', 1));
    case 'v_heq_mf2_d2'
      image = histeq(image);
      filtered_image = imdilate(medfilt2(image, [2 2]), strel('disk', 2));
    case 'v_heq_mf3_d2'
      image = histeq(image);
      filtered_image = imdilate(medfilt2(image, [3 3]), strel('disk', 2));
    case 'v_heq_d2_e2'
      image = histeq(image);
      filtered_image = imerode(imdilate(image, strel('disk', 2)), strel('disk', 2));
    case 'v_s7'
      k = ones(7);
      k = k/numel(k);
      filtered_image = conv2(image, k);
      filtered_image = filtered_image(4:end-3, 4:end-3);
    case {'v_ncut_0', 'v_ncut_0_e2', 'v_ncut_0_e4'}
      filtered_image = boundary_normalized_cuts(image, version_name, config, image_prefix);
    case 'v_neg'
      filtered_image = 1- image;
    case 'v_neg_heq_mf3'
      image_n = 1- image;
      image_n_h = histeq(image_n);
      filtered_image = medfilt2(image_n_h, [3 3]);
    case 'v_neg_heq_mf2'
      image_n = 1- image;
      image_n_h = histeq(image_n);
      filtered_image = medfilt2(image_n_h, [2 2]);
      
    otherwise
      error('Filter option not understood!!!');
  end;
  image = filtered_image;
end

if(is_verbose)
  fprintf('\nSTOP: filter_image\n');
end
return
end

