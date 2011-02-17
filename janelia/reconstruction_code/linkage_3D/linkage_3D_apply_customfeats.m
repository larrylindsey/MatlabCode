function linkage_3D_apply_customfeats (config)
% linkage_3D_gen_linkage_gph_intensity_pair_boost(config)
% generate linkage graph using boosted classifier trained on
% intensity pair histograms features
% 
% - Takes as input a stack of images and 2D superpixel and segment maps.
% - Uses trained classifier from linkage_3D_train_customfeats.m to perform 3D linkage.
% - Constructs a linkage graph, links_3D{}, with segments as nodes and 3D link
% confidences as the edge weights. This graph's link-weights are
% thresholded in proofread_linkage_superpixel.m to get the final 3D
% segmentation.
%
% AG mod based on linkage_3D_gen_linkage_gph_intensity_pair_boost
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~03202008 init code
% v1  04112008  modified for reconstruction pipeline
%

stack_config = config.stack;
seg_config = config.superpixel_2_seg;
linkage_config = config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

if( strcmp(feature_config.type, 'custom7')==0 && strcmp(feature_config.type, 'custom7_s')==0 )
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'logreg')==0)
  error('Classifier type does not match with called function. Exiting');
end;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_config.method, '/'];

junk = tree_node_w(1); %#ok<NASGU>
linkage_model = load([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
  'annotation_file', 'beta');

for iz = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(iz);

  image = get_image_from_stack(config, case_id);
  if(~isempty(stack_config.roi))
    image = image(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end;
  ims(:,:,iz) = im2double(image);
  
  segment_2D_label_map = load(sprintf([seg_dir, stack_config.image_prefix, ...
    apply_config.segmentation_suffix,'.mat'] , case_id));
  if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
    segment_2D_label_map.label_map = imresize(segment_2D_label_map.label_map, ...
      size(image), 'nearest');
  end
  segauto{iz} = segment_2D_label_map.label_map;
end

fprintf('Overlap patches .. ');tic;
ovlpList = ovlpSlice2sliceSeg(segauto);%btw-slice overlap patches
fprintf('done, time %f\n',toc);

fprintf('Calculating features .. ');tic;
[logminnorm, lightcos, minavg, litecorr] = featsLinkGrscale(ovlpList,ims);
fprintf('done - grscale, time %f    ',toc);
if 0~=strcmp(feature_config.type, 'custom7_s')
    share = featsLinkGeometry(ovlpList,segauto);
end
fprintf('done - geometry, time %f\n',toc);

% apply log regr
X = [ones(length(ovlpList),1) logminnorm' lightcos' litecorr' minavg'...
    (logminnorm.*lightcos)' (logminnorm.*litecorr)' (logminnorm.*minavg)'...
    (lightcos.*litecorr)' (lightcos.*minavg)' (litecorr.*minavg)' ];
if 0~=strcmp(feature_config.type, 'custom7_s')
    X = [X share'];
end
P = 1 ./ (1 + (exp(-X*linkage_model.beta)));
    

fprintf('Constructing linkage graph (links_3D) .. ');tic;
links_3D = [];
for iz = 1:(length(stack_config.case_ids)-1)
    zInd = [ovlpList.sliceLo]==iz;
    links_3D{iz} = [ [ovlpList(zInd).loIdx]'  [ovlpList(zInd).upIdx]'  P(zInd) ];
end
fprintf('done, time %f\n',toc);


save([get_reconstruction_dir(config), linkage_config.dir, 'links_3D', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, ...
  apply_config.segmentation_suffix, '.mat'], 'links_3D');

