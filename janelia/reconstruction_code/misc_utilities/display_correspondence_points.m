function display_correspondence_points(image_1, image_2, points_1, points_2, figure_id)
% Display correspondence points in a figure
% Input:
%   image_1, image_2      images
%   points_1, points_2    correspondence points [2xR]
%   figure_id             figure id in which to draw. If -1 then new figure
%                           is created.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  01082009    init. code
%
[height_1, width_1] = size(image_1);
scale_1 = 500/height_1;
image_1 = imresize(image_1, [500, round(scale_1*width_1)]);
[height_1, width_1] = size(image_1);
points_1 = points_1*scale_1;

[height_2, width_2] = size(image_2);
scale_2 = 500/height_2;
image_2 = imresize(image_2, [500, round(scale_2*width_2)]);
[height_2, width_2] = size(image_2);
points_2 = points_2*scale_2;

max_height = max(height_1, height_2);
display_image = zeros(max_height, width_1+width_2+10);
display_image(1:height_1, 1:width_1) = image_1;
dx = width_1 + 10;
display_image(1:height_2, dx:dx+width_2-1) = image_2;
if(figure_id<0)
  figure;
else
  figure(figure_id)
end
imshow(display_image, []);
hold on;
plot(points_1(1,:), points_1(2,:), 'ob', 'MarkerSize', 4);
plot(points_2(1,:)+dx, points_2(2,:), 'ob', 'MarkerSize', 4);
plot([points_1(1,:); points_2(1,:)+dx], [points_1(2,:); points_2(2,:)], 'r', 'LineWidth', 3);
hold off;

return
end
