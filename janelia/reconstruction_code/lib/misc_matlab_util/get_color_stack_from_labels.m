function get_color_stack_from_labels(image_tif_file, label_stack, output_tif_file)

label_stack = uint32(label_stack);
l = [];
for i = 1:size(label_stack,3)
  label_stack(:,:,i) = remove_merged_boundaries_2D(label_stack(:,:,i));
  l = [l, label_stack(:,:,i)]; %#ok<AGROW>
end
l = label2rgb(l, 'jet', [0 0 0], 'shuffle');

i=1;
fprintf('i: 1\n');
im = 1-im2double(imread(image_tif_file, i));
c = l(:,(i-1)*size(im,2)+1:i*size(im,2), :);
c = uint8(double(c).*repmat(im, [1 1 3]));
imwrite(c, output_tif_file, 'Compression', 'none');
for i = 2:size(label_stack,3)
  fprintf('i: %d\n', i);
  im = 1-im2double(imread(image_tif_file, i));
  c = l(:,(i-1)*size(im,2)+1:i*size(im,2), :);
  c = uint8(double(c).*repmat(im, [1 1 3]));
  imwrite(c, output_tif_file, 'Compression', 'none', 'WriteMode', 'append');
end
return
end
