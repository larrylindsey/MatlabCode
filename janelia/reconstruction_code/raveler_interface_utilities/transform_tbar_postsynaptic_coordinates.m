function tbar_info_str_mod = transform_tbar_postsynaptic_coordinates(...
  tbar_info_str, T, align_roi, proofreader_roi, current_z)

if(length(tbar_info_str) <= 8 || strcmp(tbar_info_str(1:5), 'T-bar')~=1)
  tbar_info_str_mod = tbar_info_str;
  return
end

tbar_info = parse_json(tbar_info_str(8:end));

tbar_str_mod = [tbar_info_str(1:7), ...
  sprintf('{"status": "%s", "partners": ', tbar_info.status)];
tbar_str_mod = [tbar_str_mod, '['];
for post_id = 1:length(tbar_info.partners)
  px = tbar_info.partners{post_id}{1}{1};
  py = tbar_info.partners{post_id}{1}{2};
  pz = tbar_info.partners{post_id}{1}{3};
  
  if(~isempty(proofreader_roi))
    py = proofreader_roi.ymax - py - proofreader_roi.ymin + 1;
  end
  
  if(~isempty(align_roi))
    py = py + align_roi.ymin-1;
    px = px + align_roi.xmin-1;
  end
  
  if(isfield(proofreader_roi, 'ymin'))
    py = py + proofreader_roi.ymin-1;
    px = px + proofreader_roi.xmin-1;
  end
  
  pxt = round(T(1,:)*[px, py, 1]');
  pyt = round(T(2,:)*[px, py, 1]');
  pzt = round(pz + current_z);
  
  tbar_str_mod = [tbar_str_mod, sprintf('[[ %d, %d, %d], %g]', pxt, pyt, pzt, ...
    tbar_info.partners{post_id}{2})]; %#ok<AGROW>
  if(post_id<length(tbar_info.partners))
    tbar_str_mod = [tbar_str_mod, ', ']; %#ok<AGROW>
  end
end
tbar_info_str_mod = [tbar_str_mod, sprintf('], "confidence": %g}', tbar_info.confidence)];

return
end
