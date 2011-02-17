function linkage_3D_train_customfeats(config)
% linkage_3D_train_intensity_pair_boost(config)
% Trains a boosted classifier to perform 3D linkage.
%
% - Features are collected for spatially overlapping 2D segments in manually
% annotated data. If two segments have the same label then their feature is
% put in the postive training set, else in the negative training set.
%
% AG mod based on linkage_3D_train_intensity_pair_boost.m
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~03202008 init code
% v1  04112008  modified for reconstruction pipeline
%

linkage_config = config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;

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

fprintf('Loading manual annotation .. ');
manual_annotation = load(train_config.manual_annotation_file);

[ny,nx,nz]=size(manual_annotation.seg);
xroi=1:nx;
yroi=1:ny;
zroi=1:nz;
if isfield(train_config,'roi')
    if isfield(train_config.roi,'xmin'), xroi=train_config.roi.xmin : train_config.roi.xmax; end
    if isfield(train_config.roi,'ymin'), yroi=train_config.roi.ymin : train_config.roi.ymax; end
    if isfield(train_config.roi,'zmin'), zroi=train_config.roi.zmin : train_config.roi.zmax; end
end

iz=1;
for iz0=zroi
    seg_gt{iz} = manual_annotation.seg(yroi,xroi,iz0);
    iz=iz+1;
end
fprintf('done\n');

fprintf('Loading images .. ');

iz=1;
for case_id = config.stack.case_ids(zroi)
  image = get_image_from_stack(config, case_id);
  if(~isempty(config.stack.roi))
    image = image(config.stack.roi.ymin:config.stack.roi.ymax, ...
      config.stack.roi.xmin:config.stack.roi.xmax);
  end;

  ims(:,:,iz) = image(yroi,xroi);
  iz=iz+1;
end
fprintf('done\n');

fprintf('Collecting linkage statistics .. ');

    %btw-slice overlap patches; use GT seg for learning
    ovlpList = ovlpSlice2sliceSeg(seg_gt);
    
    % linked or not: assign GT seg id to patch by most prevalent; assign link if same seg id
    ovlpLink01 = false(size(ovlpList));
    for io=1:length(ovlpList)
        labelsLo = seg_gt{ovlpList(io).sliceLo  }( ovlpList(io).pxls );
        labelsUp = seg_gt{ovlpList(io).sliceLo+1}( ovlpList(io).pxls ); %proofseg.seg( ovlpList(io).pxls + h*w*(ovlpList(io).sliceLo+1-1) );
        labelsLo = labelsLo( labelsLo>0 );
        labelsUp = labelsUp( labelsUp>0 );
        if isempty(labelsLo) || isempty(labelsUp) %all zeros
            ovlpLink01(io) = false;
        else
            ovlpLink01(io) = (mode(double(labelsLo)) == mode(double(labelsUp)));
        end
    end

    % feats
    [logminnorm, lightcos, minavg, litecorr] = featsLinkGrscale(ovlpList,ims);
    if 0~=strcmp(feature_config.type, 'custom7_s')
        share = featsLinkGeometry(ovlpList,seg_gt);
    end

fprintf('Training classifier .. ');
    %  train log.regr.
    Y = double(ovlpLink01)';
    X = [ones(length(ovlpList),1) logminnorm' lightcos' litecorr' minavg'...
        (logminnorm.*lightcos)' (logminnorm.*litecorr)' (logminnorm.*minavg)'...
        (lightcos.*litecorr)' (lightcos.*minavg)' (litecorr.*minavg)' ];
    if 0~=strcmp(feature_config.type, 'custom7_s')
        X = [X share'];
    end
    beta = logist2( Y, X );


annotation_file = train_config.manual_annotation_file;
save([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
  'annotation_file', 'beta');
fprintf('done\n');
