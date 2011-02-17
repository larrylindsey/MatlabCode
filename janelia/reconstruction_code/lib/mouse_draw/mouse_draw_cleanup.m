function mouse_draw_cleanup(hFig)
global mouse_draw_data

props = getappdata(hFig,'mouse_drag_parent_callbacks_0');
set(hFig,props);


if(~isempty(mouse_draw_data.hRect))
  delete(mouse_draw_data.hRect);
  mouse_draw_data.hRect = [];
end;
mouse_draw_data = [];

end