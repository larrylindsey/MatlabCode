function update_for_new_segment(new_segment_id, new_body_id, area, zpos) %#ok<INUSL>
global data
data.tmp{7}(new_segment_id, :) = [new_segment_id, area, area, zpos]
data.tmp{8}(new_segment_id, :) = 0;
data.morder(new_segment_id) = area;
data.vorder(new_segment_id) = area;
data.zmajor(new_segment_id) = zpos;
data.objk(new_segment_id,:) = [zpos, zpos];
data.major(end+1) = new_segment_id;

return
end
