function filtered_image = filter_image2(image, version_name, config, image_prefix)
% filtered_image = filter_image2(image, version_name, config, image_prefix)
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
%   filter_image
%

global config_global
if(isfield(config_global, 'DEBUG') && config_global.DEBUG)
  is_verbose = true;
else
  is_verbose = false;
end

if(is_verbose)
  fprintf('START: filter_image\n');
end

if(isempty(version_name))
  filtered_image = image;
  return;
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
  
  % For backward compatibility
  if(strcmp(version_name(1), 'v')==1)
    if(nargin==4)
      filtered_image = filter_image(image, version_name, config, image_prefix);
    elseif(nargin==3)
      filtered_image = filter_image(image, version_name, config);
    else
      filtered_image = filter_image(image, version_name);
    end
    image = filtered_image;
    continue;
  end
  
  % split the filter version name into sub-filters separated by '_'
  sub_version_names = strsplit2(version_name, '_', true);

  filtered_image = image;
  
  for s_v_id = 1:length(sub_version_names)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Don't change this part
    sub_version_name = sub_version_names{s_v_id};
    if(is_verbose)
      fprintf('sub_version_name: %s\n', sub_version_name);
    end
    image = filtered_image;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % load mitochondria detection map and set pixel within mitochondria to
    % 0
    filter_name = 'mito';
    if(apply_filter())
      mitochondria_config = config.mitochondria;
      mito = load_mitochondria_detection(...
        [get_reconstruction_dir(config), mitochondria_config.dir, ...
        image_prefix, mitochondria_config.apply.save_suffix, '.mat']);
      p = param_s();
      if(isempty(strfind(p, '_')))
        mito = mito.mitochondria_detection_confidence>param();
      else
        ss = strfind(p, '_');
        ss = ss(1);
        mito = mito.mitochondria_detection_confidence>str2double(p(1:ss-1));
        mito = filter_image2(mito, p(ss+1:end));
      end
      filtered_image = image;
      filtered_image(mito>0) = 0;
      continue;
    end
    
    
    % smooth with a square kernel of width 2*w+1
    filter_name = 's'; 
    if(apply_filter())
      w = param()*2 + 1;
      k = ones(w);
      k = k/numel(k);
      filtered_image = conv2(image, k);
      filtered_image = filtered_image(param()+1:end-param(), param()+1:end-param());
      continue;
    end
    
    % imopen
    filter_name = 'o';
    if(apply_filter())
      filtered_image = imopen(image, strel('square', param()));
      continue;
    end

    % negative
    filter_name = 'neg'; 
    if(apply_filter())
      filtered_image = 1 - image;
      continue;
    end

    % multiply
    filter_name = 'mul'; 
    if(apply_filter())
      filtered_image = image*param();
      continue;
    end

    % add
    filter_name = 'add'; 
    if(apply_filter())
      filtered_image = image + param();
      continue;
    end
    
    % div
    filter_name = 'div'; 
    if(apply_filter())
      filtered_image = image / param();
      continue;
    end
    
    %db
    filter_name = 'db';
    if(apply_filter())
      filtered_image = double(image);
      continue;
    end
    
    % nwo
    filter_name = 'nwo';
    if(apply_filter())
      filtered_image = NormalizeWithOverlap(image,param());
      continue;
    end
    
    % distr
    filter_name = 'distr';
    if(apply_filter())
      Threshedimage = zeros(size(image));
      Threshedimage(image>0.8) = 1;
      D = bwdist(Threshedimage);
      D = D./max(max(D));
      filtered_image = filtered_image-0.01*double(D);
      continue;
    end
    
    % zscore
    % must be run after double conversion
    filter_name = 'zscore';
    if(apply_filter())
      filtered_image = image - mean2(image);
      filtered_image = filtered_image./std2(filtered_image);
      continue;
    end
    
    % ps
    % must be run after double conversion
    filter_name = 'ps';
    if(apply_filter())
    p = param_s();
      psplit = strfind(p, '_');
      if(~isempty(psplit))
        Scale = str2double(p(1:psplit-1)); % 2
      else
        Scale = 2;
        psplit = 0;
      end
      filtered_image = 1./(1+exp(-((image.*Scale)-str2double(p(psplit+1:end)))));
      continue;
    end
    
    % LDA
    % TO DO: I should parameterize patch size
    filter_name = 'LDA';
    if(apply_filter())
      L = load2([get_reconstruction_dir(config), config.segmentation_2D.dir, ...
        config.LDA.dir config.LDA.train.dir 'LDAStruct' num2str(param()-1) '.mat']);
      H = reshape(L.LDAStruct.LDF,19,19);
      filtered_image = imfilter(image,H,'symmetric');     
      continue;
    end
    
    % Brute force mitochondria filter calculation
    % This filter only calculates and stores the results in T
    % but does not change filtered image.
    % Call abf filter to apply the change
    filter_name = 'cbf';
    if(apply_filter())
      FV = 7;
      L = load2([get_reconstruction_dir(config), config.segmentation_2D.dir, ...
        config.LDA.dir config.LDA.train.dir 'LDAStruct' num2str(FV) '.mat']);
      H = reshape(L.LDAStruct.LDF,19,19);
      fi = imfilter(image,H,'symmetric');
      fi = 1-fi;
      T = logical(EfficientBruteForceMitoBackComp(fi,0.9,param()));
      TT = logical(EfficientBruteForceMitoBackComp(fi,0.8,param()));
      CBFtmpFilter = T | TT;
      continue;
    end
    
    % Brute force mitochondria filter application
    % Apply before neg filter
    filter_name = 'abf';
    if(apply_filter())
      filtered_image = image;
      filtered_image(CBFtmpFilter) = 0;
      continue;
    end
    
    % rescale
    filter_name = 'rescale';
    if(apply_filter())
      filtered_image = image  - min(min(image));
      filtered_image = filtered_image./max(max(filtered_image));
      continue;
    end
 
      
    % median filtering
    filter_name = 'mf';
    if(apply_filter())
      filtered_image = medfilt2(image, param()*[1 1], 'symmetric');
      continue;
    end

    % adaptive histogram equalization
    filter_name = 'aheq'; 
    if(apply_filter())
      filtered_image = adapthisteq(image);
      continue;
    end
    
    % histogram equalization - excluding outlier intensities possibly
    % belonging to folds and tears
    filter_name = 'heq2'; 
    if(apply_filter())
      filtered_image = histeq2(image(image>0.1 & image<0.9), image);
      continue;
    end
    % histogram equalization
    filter_name = 'heq'; 
    if(apply_filter())
      filtered_image = histeq(image);
      continue;
    end
    
    % convert to gray
    filter_name = 'g'; 
    if(apply_filter())
      if(size(image,3)==3)
        filtered_image = rgb2gray(image);
      else
        filtered_image = image;
      end
      continue;
    end
    
    % erode only zero-valued pixels
    filter_name = 'ez'; 
    if(apply_filter())
      temp_image = image==0;
      temp_image = imdilate(temp_image, strel('disk', param()));
      filtered_image = image;
      filtered_image(temp_image>0) = 0;
      continue;
    end
    
    % erosion
    filter_name = 'e'; 
    if(apply_filter())
      filtered_image = imerode(image, strel('disk', param()));
      continue;
    end

    % dilate only onto zero-valued pixels
    filter_name = 'dz'; 
    if(apply_filter())
      temp_image = imdilate(image, strel('disk', param()));
      filtered_image = image;
      filtered_image(image==0) = temp_image(image==0);
      continue;
    end
    
    % dilation
    filter_name = 'd'; 
    if(apply_filter())
      filtered_image = imdilate(image, strel('disk', param()));
      continue;
    end

    % imclose
    filter_name = 'c';
    if(apply_filter())
      filtered_image = imclose(image, strel('square', param()));
      continue;
    end

    % bwareaopen
    filter_name = 'a';
    if(apply_filter())
      filtered_image = bwareaopen(image, param());
      continue;
    end

    % threshold
    filter_name = 't';
    if(apply_filter())
      filtered_image = double(image>param());
      continue;
    end
    
    % lower bound
    filter_name = 'lb';
    if(apply_filter())
      filtered_image = image;
      filtered_image(image<param()) = param();
      continue;
    end
    
    % upper bound
    filter_name = 'ub';
    if(apply_filter())
      filtered_image = image;
      filtered_image(image>param()) = param();
      continue;
    end
    
    % load the BEL boundary map
    filter_name = 'BEL';
    if(apply_filter())
      stack_config = config.stack;
      seg_config = config.segmentation_2D;
      BEL_dir = [get_reconstruction_dir(config), seg_config.dir, config.BEL.dir];
      BEL_name = [BEL_dir, image_prefix, '.BEL.', param_s(), '.mat'];
      BEL = load2(BEL_name, 'boundary');
      if(~isempty(stack_config.roi))
        boundary = BEL.boundary(stack_config.roi.ymin:stack_config.roi.ymax, ...
          stack_config.roi.xmin:stack_config.roi.xmax);
      else
        boundary = BEL.boundary;
      end
      if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
        filtered_image = imresize(boundary, 1/stack_config.segmentation_scale, ...
          stack_config.segmentation_scaling_method);
      else
        filtered_image = boundary;
      end
      continue;
    end
    
    error('filter not recognized.');
  end
  
  image = filtered_image;
end
  function apply = apply_filter()
    apply = strncmp(sub_version_name, filter_name, length(filter_name))==1;
    return
  end

  function p = param()
    p = str2double(sub_version_name(length(filter_name)+1:end));
    return
  end

  function p = param_s()
    p = sub_version_name(length(filter_name)+1:end);
    p = regexprep(p, '(?<char>[^\\])\', '$<char>');
    return
  end

if(is_verbose)
  fprintf('STOP: filter_image\n');
end
return
end

