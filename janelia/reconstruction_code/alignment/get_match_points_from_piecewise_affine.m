function matches = get_match_points_from_piecewise_affine(...
  transforms_tp, fold_mask_1, fold_mask_2, min_n_match)
% matches = get_match_points_from_piecewise_affine(...
%   transforms_tp, fold_mask_1, fold_mask_2, min_n_match)
% Given a piecewise affine transformation transforms_tp from image1 to
% image2, get a set of correspondences from image1 to image2.
% If fold_mask_1 and fold_mask_2 to are specified then the correspondences
% are partitioned accordingly.
% If min_n_match is specified then atleast these many correspondence pairs
% are returned.
% Input:
%   transforms_tp     piecewise affine transform structure
%     .map_mask       map of R transform domains (MxN uint16 matrix)
%     .transform      affine transforms (6xR double matrix)
%   fold_mask_1       fold mask of connected regions in image1 (MxN uint8
%                       matrix) 
%   fold_mask_2       fold mask of connected regions in image2 (M2xN2 uint8   
%                       matrix) 
%   min_n_match       minimum number of correspondence pairs in image1
%                       desired (scalar double)
%
% Output:
%   matches           struct of correspondences
%     .points_1       [x;y] points in image1
%     .points_2       [x;y] points in image2
%     .set_ids_1      partition id points_1 if fold_mask_1 and fold_mask_2
%                       are specified
%     .set_ids_2      partition id points_2 if fold_mask_1 and fold_mask_2
%                       are specified
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0    01082009    init. code
%

global config_global
TRANSFORMATION_ID_OFFSET = config_global.TRANSFORMATION_ID_OFFSET;

if(nargin>=2)
  [height_2, width_2] = size(fold_mask_2);
  is_specified_fold_mask = true;
else
  height_2 = inf;
  width_2 = inf;
  is_specified_fold_mask = false;
end

if(nargin<4)
  min_n_match = 0;
end

% Get a correspondence pair for each region in the piecewise transform.
matches.points_1 = [];
matches.points_2 = [];
if(is_specified_fold_mask)
  matches.set_ids_1 = [];
  matches.set_ids_2 = [];
end
n_transform_txt = size(transforms_tp.transforms,2);
print_debug('no. of transforms in transform struct = %d\n', n_transform_txt);
n_transform_mask = max(transforms_tp.map_mask(:)) - TRANSFORMATION_ID_OFFSET + 1;
print_debug('no. of transforms in map_mask patch ids = %d\n', n_transform_mask);
n_transform = min(n_transform_mask, n_transform_txt);
for transform_id = 1:n_transform
  transform_mask = transforms_tp.map_mask == (transform_id+TRANSFORMATION_ID_OFFSET-1);
  [py, px] = find(transform_mask);
  x1 = round(mean(px));
  y1 = round(mean(py));
  transform = reshape(transforms_tp.transforms(:,transform_id), [2 3]);
  x2 = round(transform(1,:)*[x1; y1; 1]);
  y2 = round(transform(2,:)*[x1; y1; 1]);
  if(x2<=0 || x2>width_2 || y2<=0 || y2>height_2)
    continue;
  end
  matches.points_1(:,end+1) = [x1; y1];
  matches.points_2(:,end+1) = [x2; y2];
  if(is_specified_fold_mask)
    matches.set_ids_1(end+1) = fold_mask_1(y1, x1);
    matches.set_ids_2(end+1) = fold_mask_2(y2, x2);
  end
end

% Make sure we have minimum number of matches
if(~is_specified_fold_mask)
  n_match_to_get = min_n_match - size(matches.points_1, 2);
  if(n_match_to_get<=0)
    return;
  end
  mask = imerode(transforms_tp.map_mask>0, strel('disk', 10));
  [py, px] = find(mask);
  if(isempty(px))
    [py, px] = find(transforms_tp.map_mask);
  end
  random_point_ids = 1 + round((length(px)-1)*rand([n_match_to_get,1]));
  for i = random_point_ids
    x1 = px(i);
    y1 = py(i);
    transform = reshape(transforms_tp.transforms(:,transform_id), [2 3]);
    x2 = round(transform(1,:)*[x1; y1; 1]);
    y2 = round(transform(2,:)*[x1; y1; 1]);
    if(x2<=0 || x2>width_2 || y2<=0 || y2>height_2)
      continue;
    end
    matches.points_1(:,end+1) = [x1; y1];
    matches.points_2(:,end+1) = [x2; y2];
  end
end

n_set_1 = max(fold_mask_1(:));
for set_id = 1:n_set_1
  n_match_to_get = min_n_match - sum(matches.set_ids_1==set_id);
  if(n_match_to_get<=0)
    continue;
  end
  mask = imerode((transforms_tp.map_mask>0) & (fold_mask_1==set_id), ...
    strel('disk', 10));
  [py, px] = find(mask);
  if(isempty(px))
    [py, px] = find((transforms_tp.map_mask>0) & (fold_mask_1==set_id));
  end
  if(isempty(px))
    continue;
  end
  random_point_ids = 1 + round((length(px)-1)*rand([n_match_to_get,1]));
  for i = 1:length(random_point_ids)
    x1 = px(random_point_ids(i));
    y1 = py(random_point_ids(i));
    transform = reshape(transforms_tp.transforms(:,transform_id), [2 3]);
    x2 = round(transform(1,:)*[x1; y1; 1]);
    y2 = round(transform(2,:)*[x1; y1; 1]);
    if(x2<=0 || x2>width_2 || y2<=0 || y2>height_2)
      continue;
    end
    matches.points_1(:,end+1) = [x1; y1];
    matches.points_2(:,end+1) = [x2; y2];
    matches.set_ids_1(end+1) = fold_mask_1(y1, x1);
    matches.set_ids_2(end+1) = fold_mask_2(y2, x2);
  end
end

return
end
