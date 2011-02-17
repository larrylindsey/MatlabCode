function image_filter(image_in_file_name, image_out_file_name, filter_config, is_verbose)
% image_filter(image_in_file_name, image_out_file_name, filter_config, is_verbose)
% Read in an image image_in_file_name, filter it based on filter config and
% then write it out to image_out_file_name
%
% is_verbose    print messages
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

if(nargin==1 && strcmp(image_in_file_name, '-h')==1)
  fprintf('image_filter image_in_file_name image_out_file_name filter_config [is_verbose]\n');
  fprintf('  is_verbose    print messages [true]\n');
  fprintf('  filter_config:\n');
  fprintf('    s<w>      smooth with a square kernel of width 2*w+1\n');
  fprintf('    o<w>      morphological open with a square kernel w\n');
  fprintf('    neg<w>    w-I\n');
  fprintf('    add<w>    w+I\n');
  fprintf('    mul<w>    w*I\n');
  fprintf('    div<w>    I/w\n');
  fprintf('    rcp<w>    w/I\n');
  fprintf('    mf<w>     median filter with kernel of width w\n');
  fprintf('    aheq      adaptive histogram equalization\n');
  fprintf('    heq2      histogram equalization while ignoring intensities less than 0.1 and greater than 0.9\n');
  fprintf('    heq       histogram equalization\n');
  fprintf('    g         convert to grayscale\n');
  fprintf('    cd        convert to double - no normalization\n');
  fprintf('    cdn       convert to double and normalize to range [0,1]\n');
  fprintf('    cuc       convert to unsigned char\n');
  fprintf('    ez<w>     morphological erosion with disk of width w only 0 valued pixels\n');
  fprintf('    e<w>      morphological erosion with disk of width w\n');
  fprintf('    dz<w>     morphological dilation with disk of width w only onto 0 valued pixels\n');
  fprintf('    e<w>      morphological dilation with disk of width w\n');
  fprintf('    c<w>      morphological close with square of width w\n');
  fprintf('    a<w>      remove connected components of area less than w\n');
  fprintf('    t<w>      I>w?1:0\n');
  return;
end

if(isdeployed)
  if(nargin<=3)
    is_verbose = true;
  else
    is_verbose = eval(is_verbose);
  end
else
  if(nargin<=3)
    is_verbose = true;
  end
end

if(is_verbose)
  fprintf('START: filter_image\n');
end

image = imread(image_in_file_name);

if(isempty(filter_config))
  imwrite(image, image_out_file_name);
  return;
end

% split the filter version name into sub-filters separated by '_'
sub_filter_configs = strsplit2(filter_config, '_', true);

filtered_image = image;
for s_v_id = 1:length(sub_filter_configs)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Don't change this part
  fc = sub_filter_configs{s_v_id};
  if(is_verbose)
    fprintf('fc: %s\n', fc);
  end
  image = filtered_image;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % smooth with a square kernel of width 2*w+1
  fn = 's';
  if(apply_filter(fc,fn))
    w = param(fc,fn)*2 + 1;
    k = ones(w);
    k = k/numel(k);
    filtered_image = conv2(image, k);
    filtered_image = filtered_image(param(fc,fn)+1:end-param(fc,fn), ...
      param(fc,fn)+1:end-param(fc,fn));
    continue;
  end
  
  % imopen
  fn = 'o';
  if(apply_filter(fc,fn))
    filtered_image = imopen(image, strel('square', param(fc,fn)));
    continue;
  end
  
  % negative
  fn = 'neg';
  if(apply_filter(fc,fn))
    filtered_image = param(fc,fn) - image;
    continue;
  end
  
  % add
  fn = 'add';
  if(apply_filter(fc,fn))
    filtered_image = param(fc,fn) + image;
    continue;
  end
  
  % multiply
  fn = 'mul';
  if(apply_filter(fc,fn))
    filtered_image = image*param(fc,fn);
    continue;
  end
  
  % divide
  fn = 'div';
  if(apply_filter(fc,fn))
    filtered_image = image/param(fc,fn);
    continue;
  end
  
  % divide
  fn = 'rcp';
  if(apply_filter(fc,fn))
    filtered_image = param(fc,fn)./image;
    continue;
  end
  
  % median filtering
  fn = 'mf';
  if(apply_filter(fc,fn))
    filtered_image = medfilt2(image, param(fc,fn)*[1 1], 'symmetric');
    continue;
  end
  
  % adaptive histogram equalization
  fn = 'aheq';
  if(apply_filter(fc,fn))
    filtered_image = adapthisteq(image);
    continue;
  end
  
  % histogram equalization - excluding outlier intensities possibly
  % belonging to folds and tears
  fn = 'heq2';
  if(apply_filter(fc,fn))
    filtered_image = histeq2(image(image>0.1 & image<0.9), image);
    continue;
  end
  % histogram equalization
  fn = 'heq';
  if(apply_filter(fc,fn))
    filtered_image = histeq(image);
    continue;
  end
  
  % convert to gray
  fn = 'g';
  if(apply_filter(fc,fn))
    if(size(image,3)==3)
      filtered_image = rgb2gray(image);
    else
      filtered_image = image;
    end
    continue;
  end

  % convert to double - no normalization
  fn = 'cd';
  if(apply_filter(fc,fn))
    filtered_image = double(image);
    continue;
  end

  % convert to unsigned char
  fn = 'cuc';
  if(apply_filter(fc,fn))
    filtered_image = uint8(image);
    continue;
  end
  
  % erode only zero-valued pixels
  fn = 'ez';
  if(apply_filter(fc,fn))
    temp_image = image==0;
    temp_image = imdilate(temp_image, strel('disk', param(fc,fn)));
    filtered_image = image;
    filtered_image(temp_image>0) = 0;
    continue;
  end
  
  % erosion
  fn = 'e';
  if(apply_filter(fc,fn))
    filtered_image = imerode(image, strel('disk', param(fc,fn)));
    continue;
  end
  
  % dilate only onto zero-valued pixels
  fn = 'dz';
  if(apply_filter(fc,fn))
    temp_image = imdilate(image, strel('disk', param(fc,fn)));
    filtered_image = image;
    filtered_image(image==0) = temp_image(image==0);
    continue;
  end
  
  % dilation
  fn = 'd';
  if(apply_filter(fc,fn))
    filtered_image = imdilate(image, strel('disk', param(fc,fn)));
    continue;
  end
  
  % imclose
  fn = 'c';
  if(apply_filter(fc,fn))
    filtered_image = imclose(image, strel('square', param(fc,fn)));
    continue;
  end
  
  % bwareaopen
  fn = 'a';
  if(apply_filter(fc,fn))
    filtered_image = bwareaopen(image, param(fc,fn));
    continue;
  end
  
  % threshold
  fn = 't';
  if(apply_filter(fc,fn))
    filtered_image = image>param(fc,fn);
    continue;
  end
end

imwrite(filtered_image, image_out_file_name);

if(is_verbose)
  fprintf('STOP: filter_image\n');
end
return
end


function apply = apply_filter(fc,fn)
apply = strncmp(fc,fn, length(fn))==1;
return
end

function p = param(fc,fn)
p = str2double(fc(length(fn)+1:end));
return
end
