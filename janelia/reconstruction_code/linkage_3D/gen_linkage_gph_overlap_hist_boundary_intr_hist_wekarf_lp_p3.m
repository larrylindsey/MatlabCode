function gen_linkage_gph_overlap_hist_boundary_intr_hist_wekarf_lp_p3(config)
% gen_linkage_gph_overlap_hist_boundary_intr_hist_lp_p3(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  12122008  init code
%

if(isfield(config.align.linkage_align, 'method') && ...
    ~isempty(config.align.linkage_align.method))
  linkage_align_method = config.align.linkage_align.method;
else
  linkage_align_method = 'prealign';
end
switch(linkage_align_method)
  case {'prealign', 'SIFT'}
    error('To Do');
%     gen_link_gph_overlap_hist_boundary_intr_hist_lp_prealign_sift(config)
  case 'deformable_mesh'
    gen_link_gph_overlap_hist_boundary_intr_hist_lp_wekarf_dmesh(config)
  otherwise
    error('Linkage section alignment not recognized');
end
return
end
