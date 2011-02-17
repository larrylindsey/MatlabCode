function body_stack = get_3D_stack_from_ws_seg_raw(ws_file_name, seg_file_name)

ws = fread_raw_array_mex(ws_file_name);
seg = fread_raw_array_mex(seg_file_name);
size_ws = size(ws);

body_stack = uint32(zeros(size_ws(2), size_ws(1), size_ws(3)));
for z = 1:size_ws(3)
    body_stack(:,:,z) = seg(1+ws(:,:,z)');
end
return;
end
