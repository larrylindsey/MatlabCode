function image = plot_markers(image, px, py, marker)
% image = plot_markers(image, px, py, marker)
% plot markers on images
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

switch(marker.color)
  case 'r'
    color = [1 0 0];
  case 'g'
    color = [0 1 0];
  case 'b'
    color = [0 0 1];
  case 'c'
    color = [0 1 1];
  case 'm'
    color = [1 0 1];
  case 'y'
    color = [1 1 0];
  case 'k'
    color = [0 0 0];
  case 'w'
    color = [1 1 1];
  otherwise
    error('Color not recognized');
end

switch(marker.type)
  case '.'
    px1 = max(px-1,1); px2 = min(px+1,size(image,2));
    py1 = max(py-1,1); py2 = min(py+1,size(image,1));
    p = [px, py; px1, py; px2, py; px, py1; px, py2; ...
      px1,py1; px1,py2; px2,py1; px2,py2];
  otherwise
    error('Marker type not recognized');
end

id = (p(:,1)-1)*size(image,1) + p(:,2);
for i = 1:3
  temp = image(:,:,i);
  temp(id) = color(i);
  image(:,:,i)=temp;
end;

return
end
