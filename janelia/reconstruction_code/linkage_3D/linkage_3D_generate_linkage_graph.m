function linkage_3D_generate_linkage_graph(config)
% linkage_3D_generate_linkage_graph(config)
% Calls a function to generate 3D linkage graph based on the value of
% config.linkage.feature.type and config.linkage.model.type
%
% Pratim Ghosh
% Summer Intern
% University of California, Santa Barbara
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  06162008  init code
%

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
        linkage_3D_gen_linkage_gph_intensity_pair_boost(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'intensity_pair_hist_v2', 'intensity_pair_hist_v2b', ...
      'intensity_pair_hist_v2c', 'intensity_pair_hist_v2d'}
    switch(config.linkage.model.type)
      case 'boost'
        linkage_3D_gen_linkage_gph_intensity_pair_v2_boost(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'custom7', 'custom7_s'}
    switch(config.linkage.model.type)
      case 'logreg'
        linkage_3D_apply_customfeats(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'overlap_area'
    switch(config.linkage.model.type)
      case 'boost'
        linkage_3D_gen_linkage_gph_overlap_area_boost(config);
      case 'boost_lp_p3'
        linkage_3D_gen_linkage_gph_overlap_area_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'normalize_overlap_area'
    switch(config.linkage.model.type)
      case 'boost'
        switch(linkage_align_method)
          case {'prealign', 'SIFT'}
            gen_linkage_gph_normalize_overlap_area_boost_prealign_sift(config)
          otherwise
            error('Linkage section alignment not recognized');
        end
      otherwise
        error('Linkage model type not recognized');
    end
  case {'overlap_hist', 'overlap_hist_shrink_LS'}
    switch(config.linkage.model.type)
      case 'boost'
        gen_linkage_gph_overlap_hist_boost(config);
      otherwise
        error('Linkage model type not recognized');
    end
    
  case 'symm_diff'
    switch(config.linkage.model.type)
      case 'lp_p3'
        linkage_3D_gen_linkage_gph_symm_diff_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'symm_diff_mean_boundary'
    switch(config.linkage.model.type)
      case 'lp_p3'
        linkage_3D_gen_linkage_gph_symm_diff_mean_boundary_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'overlap_area_boundary_hist'
    switch(config.linkage.model.type)
      case 'boost_lp_p3'
        linkage_3D_gen_linkage_gph_overlap_area_boundary_hist_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case 'overlap_area_boundary_interior_hist'
    switch(config.linkage.model.type)
      case 'boost_lp_p3'
        gen_linkage_gph_overlap_area_boundary_intr_hist_lp_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'overlap_hist_boundary_interior_hist', ...
      'overlap_hist_boundary_interior_hist_shrink_LS'}
    switch(config.linkage.model.type)
      case 'wekarf_lp_p3'
        gen_linkage_gph_overlap_hist_boundary_intr_hist_wekarf_lp_p3(config);
      case 'boost_lp_p3'
        gen_linkage_gph_overlap_hist_boundary_intr_hist_lp_p3(config);
      case 'boost_sdp'
        gen_linkage_gph_overlap_hist_boundary_intr_hist_sdp_p3(config);
      case 'boost_ncuts'
        gen_linkage_gph_overlap_hist_boundary_intr_hist_ncut_p3(config);
      otherwise
        error('Linkage model type not recognized');
    end
    
  otherwise
    error('Linkage feature type not recognized');
end

end