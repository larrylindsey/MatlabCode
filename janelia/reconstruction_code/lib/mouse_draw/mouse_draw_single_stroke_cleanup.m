function mouse_draw_single_stroke_cleanup()
global mouse_draw_data

if(exist('mouse_draw_data', 'var') && ~isempty(mouse_draw_data) && ~isempty(mouse_draw_data.hRect))
  try
    delete(mouse_draw_data.hRect);
  catch
    do_nothing = true;
  end;
  mouse_draw_data.hRect = [];
end;

mouse_draw_data = [];

end