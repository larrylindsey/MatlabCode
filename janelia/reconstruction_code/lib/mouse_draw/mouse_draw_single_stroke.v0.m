function mouse_draw_single_stroke(fig, hAxis, width, height, brush_size, brush_color)
% mouse_draw(fig, hAxis, width, height, brush_size, brush_color)
% draw with mouse like a paint brush - capture only one continuous stroke.
% Different from mouse_draw() which captures multiple strokes until stopped
% 
% inputs:
%   fig           handle to the figure
%   hAxis         handle to the current axis
%   width, height size of the image to be painted on (in current axis)
%   brush_size    rectangular brush of size (2*brush_size)x(2*brush_size)
%   brush_color   color to be used when drawing on figure
%
% output
%   mouse_draw_data - a global variable
%     hRect             handles to the drawn rectangles
%     mouse_points      Nx2 matrix of <x,y> coords.
%
% Call mouse_draw_cleanup() after saving the mouse_points for housekeeping.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
% 04022008
%

if(~exist('brush_color', 'var'))
  brush_color = 'b';
end;

global mouse_draw_data

mouse_draw_data.brush_size = brush_size;
% mouse_draw_data.hRect = [];
% mouse_draw_data.mouse_points = [];
mouse_draw_data.brush_color = brush_color;
axis_position = get(hAxis, 'Position');
mouse_draw_data.x_min = axis_position(1);
mouse_draw_data.y_min = axis_position(2);
mouse_draw_data.x_max = mouse_draw_data.x_min + axis_position(3);
mouse_draw_data.y_max = mouse_draw_data.y_min + axis_position(4);
mouse_draw_data.x_factor = axis_position(3);
mouse_draw_data.y_factor = axis_position(4);
mouse_draw_data.width = width;
mouse_draw_data.height = height;

set(fig, 'Pointer', 'crosshair');

% get the values and store them in the figure's appdata
props.WindowButtonMotionFcn = get(fig,'WindowButtonMotionFcn');
props.WindowButtonUpFcn = get(fig,'WindowButtonUpFcn');
setappdata(fig,'mouse_drag_parent_callbacks_1',props);

% set the new values for the WindowButtonMotionFcn and
% WindowButtonUpFcn
set(fig,'WindowButtonMotionFcn',{@wbm})
set(fig,'WindowButtonUpFcn',{@wbu})

global mouse_draw_data

% get current point and transform to image axis
c_p = get(fig, 'CurrentPoint');
if(c_p(1)<mouse_draw_data.x_min || c_p(1)>mouse_draw_data.x_max)
  return;
end;
if(c_p(2)<mouse_draw_data.y_min || c_p(2)>mouse_draw_data.y_max)
  return;
end;
mouse_draw_data.mouse_points(end+1,1) = mouse_draw_data.width*(c_p(1)-mouse_draw_data.x_min)/mouse_draw_data.x_factor;
mouse_draw_data.mouse_points(end,2) = mouse_draw_data.height*(1-(c_p(2)-mouse_draw_data.y_min)/mouse_draw_data.y_factor);

% draw rectangle
hold on;
mouse_draw_data.hRect(end+1) = rectangle('Position', [mouse_draw_data.mouse_points(end,1)-mouse_draw_data.brush_size, ...
  mouse_draw_data.mouse_points(end,2)-mouse_draw_data.brush_size, 2*mouse_draw_data.brush_size, ...
  2*mouse_draw_data.brush_size], 'FaceColor', mouse_draw_data.brush_color);
hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function called upon mouse movement during the button press
% Accumulate mouse locations into mouse_draw_data.mouse_points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbm(h,evd)
% disp('motion')
global mouse_draw_data

% get current point and transform to image axis
c_p = get(h, 'CurrentPoint');
if(c_p(1)<mouse_draw_data.x_min || c_p(1)>mouse_draw_data.x_max)
  return;
end;
if(c_p(2)<mouse_draw_data.y_min || c_p(2)>mouse_draw_data.y_max)
  return;
end;
mouse_draw_data.mouse_points(end+1,1) = mouse_draw_data.width*(c_p(1)-mouse_draw_data.x_min)/mouse_draw_data.x_factor;
mouse_draw_data.mouse_points(end,2) = mouse_draw_data.height*(1-(c_p(2)-mouse_draw_data.y_min)/mouse_draw_data.y_factor);

% draw rectangle
hold on;
mouse_draw_data.hRect(end+1) = rectangle('Position', [mouse_draw_data.mouse_points(end,1)-mouse_draw_data.brush_size, ...
  mouse_draw_data.mouse_points(end,2)-mouse_draw_data.brush_size, 2*mouse_draw_data.brush_size, ...
  2*mouse_draw_data.brush_size], 'FaceColor', mouse_draw_data.brush_color);
hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called upon release of mouse button - removes the button motion and up
% hooks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbu(h,evd)
% executes when the mouse button is released
% disp('up')

% get the properties and restore them
props = getappdata(h,'mouse_drag_parent_callbacks_1');
set(h,props);
end
