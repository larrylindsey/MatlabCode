function plot_segment_boundaries(image, label_map, scale, figure_id)
% figure;
% imshow(image);
% [py,px]=find(label_map==0);
% hold on;
% plot(px, py, '.');
% hold off;

if(nargin<3)
  scale = 1;
end

if(nargin<4)
  figure;
else
  figure(figure_id);
end
if(scale==1)
  imshow(image);
else
  imshow(imresize(image, 1/scale));
end
[py,px]=find(label_map==0);
hold on;
plot(px/scale, py/scale, '.');
hold off;

return
end
