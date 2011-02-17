function linkage_3D_train(config, sub_module_id)
% linkage_3D_train(config)
% Calls a function to train 3D linkage algorithm based on the value of
% config.linkage.feature.type and config.linkage.model.type
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  06162008  init code
%

if(config.is_verbose)
  fprintf('START: linkage_3D_train\n');
end

if(isfield(config.align.linkage_align, 'method') && ...
    ~isempty(config.align.linkage_align.method))
  linkage_align_method = config.align.linkage_align.method;
else
  linkage_align_method = 'prealign';
end

switch(config.linkage.feature.type)
  case {'intensity_pair_hist'}
    switch(config.linkage.model.type)
      case 'boost'
        linkage_3D_train_intensity_pair_boost(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'intensity_pair_hist_v2', 'intensity_pair_hist_v2b', ...
      'intensity_pair_hist_v2c', 'intensity_pair_hist_v2d'}
    switch(config.linkage.model.type)
      case 'boost'
        linkage_3D_train_intensity_pair_v2_boost(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'custom7', 'custom7_s'}
    switch(config.linkage.model.type)
      case 'logreg'
        linkage_3D_train_customfeats(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'overlap_area'
    switch(config.linkage.model.type)
      case 'boost'
        linkage_3D_train_overlap_area_boost(config);
      case 'boost_lp_p3'
        linkage_3D_train_overlap_area_boost_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'normalize_overlap_area'
    switch(config.linkage.model.type)
      case 'boost'
        linkage_3D_train_normalize_overlap_area_boost(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'overlap_hist', 'overlap_hist_shrink_LS'}
    switch(sub_module_id)
      case 10
        switch(linkage_align_method)
          case {'prealign', 'SIFT'}
            error('No need to module 710.1 for pre-aligned data. Run module 750.5 directly');
          case 'deformable_mesh'
            linkage_train_collect_feat_overlap_hist_multi_tile(config);
          otherwise
            error('Linkage section alignment not recognized');
        end
      case 50
        switch(linkage_align_method)
          case {'prealign', 'SIFT'}
            linkage_3D_train_overlap_hist_boost_prealign(config);
          case 'deformable_mesh'
            linkage_train_overlap_hist_boost_multi_tile(config);
          otherwise
            error('Linkage section alignment not recognized');
        end
    end
    
  case 'symm_diff'
    switch(config.linkage.model.type)
      case 'lp_p3'
        linkage_3D_train_symm_diff_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'symm_diff_mean_boundary'
    switch(config.linkage.model.type)
      case 'lp_p3'
        linkage_3D_train_symm_diff_mean_boundary_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'overlap_area_boundary_hist'
    switch(config.linkage.model.type)
      case 'boost_lp_p3'
        linkage_3D_train_overlap_area_boundary_hist_boost_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'overlap_area_boundary_interior_hist'
    switch(config.linkage.model.type)
      case 'boost_lp_p3'
        linkage_3D_train_overlap_area_boundary_intr_hist_boost_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'overlap_hist_boundary_interior_hist', ...
      'overlap_hist_boundary_interior_hist_shrink_LS'}
    switch(sub_module_id)
      case 10
        switch(linkage_align_method)
          case {'prealign', 'SIFT'}
            error('No need to module 710.1 for pre-aligned data. Run module 750.5 directly');
          case 'deformable_mesh'
            linkage_train_collect_feat_overlap_hist_lp_multi_tile(config);
          otherwise
            error('Linkage section alignment not recognized');
        end
      case 50
        switch(linkage_align_method)
          case {'prealign', 'SIFT'}
            linkage_3D_train_overlap_hist_boundary_intr_hist_boost_lp_p3(config);
          case 'deformable_mesh'
            linkage_train_overlap_hist_boost_lp_multi_tile(config);
          otherwise
            error('Linkage section alignment not recognized');
        end
    end
%     switch(config.linkage.model.type)
%       case {'boost_lp_p3', 'boost_sdp'}
%         linkage_3D_train_overlap_hist_boundary_intr_hist_boost_lp_p3(config);
%       case 'boost_ncuts'
%         linkage_3D_train_overlap_hist_boundary_intr_hist_boost_ncuts(config);
%       otherwise
%         error('Linkage model type not recognized');
%     end
    
  otherwise
    error('Linkage feature type not recognized');
end

if(config.is_verbose)
  fprintf('STOP: linkage_3D_train\n');
end
return
end