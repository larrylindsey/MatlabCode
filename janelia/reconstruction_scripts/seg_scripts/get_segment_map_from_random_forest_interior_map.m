function [label_map, watershed_label_map] = ...
  get_segment_map_from_random_forest_interior_map(...
  image, interior_map, varargin)
% label_map = get_segment_map_from_random_forest_boundary_map(image, interior_map)
% Compute a segment label map given an image and its boundary map detected
% using Random Forest classifier.
%
% Inputs:
%   image           grayscale, MxN uint8
%   interior_map    MxN uint8, interior points are labelled 255 and
%                     boundary points are labelled 0.
% Optional Inputs:
%   is_verbose                whether to print intermediate messages [true].
%   is_verbose_figures        whether to show intermediate results as figures
%                             [true].
%   display_scale             scale down images for display [1].
%   mitochondria_map          mitochondria map MxN, mitochondria pixels
%                             marked to >0, everything else set to 0.
%   interior_area_threshold   islands in interior map below this threshold
%                             are removed before seeding [10].
% Output:
%   label_map       MxN uint32 segment label map with 1-pixel thin
%                     boundaries between segments labelled with 0's.
%   watershed_label_map  Watershed segment map
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

is_verbose = true;
is_verbose_figures = true;
display_scale = 1;
mitochondria_map = [];
interior_area_threshold = 0;
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'is_verbose'
      is_verbose = varargin{i+1};
    case 'is_verbose_figures'
      is_verbose_figures = varargin{i+1};
    case 'display_scale'
      display_scale = varargin{i+1};
    case 'mitochondria_map'
      mitochondria_map = varargin{i+1};
    case 'interior_area_threshold'
      interior_area_threshold = varargin{i+1};
  end
end

if(display_scale<=0)
  display_scale = 1;
end

figure(1);
if(is_verbose_figures)
%   imshow(imresize(image, 1/display_scale));
%   title('image');

  figure(2);
  imshow(imresize(interior_map, 1/display_scale, 'nearest'));
  title('interior map');
end

if(is_verbose)
  fprintf('Computing initial boundary map\n');
end
boundary_map = filter_image2(image, 'db_div255_mf3_neg');
if(is_verbose_figures)
%   figure(3);
%   imshow(imresize(boundary_map, 1/display_scale));
%   title('boundary map');
end

if(is_verbose)
  fprintf('Removing small islands (area<%d) from interior_map\n', ...
    interior_area_threshold);
end
interior_map = bwareaopen(interior_map, interior_area_threshold);
if(is_verbose_figures)
  figure(4);
  imshow(imresize(interior_map, 1/display_scale, 'nearest'));
  title('modified interior map');
end

if(is_verbose)
  fprintf('Computing seeded boundary map using interior_map\n');
end
seeded_boundary_map = imimposemin(boundary_map, interior_map);
if(is_verbose_figures)
%   figure(5);
%   imshow(imresize(seeded_boundary_map, 1/display_scale));
%   title('seeded boundary map');
end

if(is_verbose)
  fprintf('Computing watershed segment map\n');
end
watershed_label_map = watershed(seeded_boundary_map);
if(is_verbose_figures)
  plot_segment_boundaries(image, watershed_label_map, display_scale, 6);
  title('watershed map');
end

if(~isempty(mitochondria_map))
  if(is_verbose)
    fprintf('Computing high-confidence mitochondria regions\n');
  end
  mitochondria_map = double(mitochondria_map>0);
  mitochondria_map = imerode(mitochondria_map, strel('disk', 10));
  mitochondria_map = bwareaopen(mitochondria_map, 40);
%   mitochondria_map = imerode(mitochondria_map, strel('disk', 20));
  boundary_map_smito = boundary_map;
  boundary_map_smito(mitochondria_map>0) = 0;
  label_map = compute_segmentation_from_superpixels_with_min_area_c(...
    watershed_label_map, boundary_map_smito, 0.01, 300); 
end

if(is_verbose_figures && nnz(label_map~=watershed_label_map)>0)
  plot_segment_boundaries(image, label_map, display_scale, 7);
  title('final label map');
end
label_map = uint32(label_map);
watershed_label_map = uint32(watershed_label_map);
return
end
