function im = setImBorder(im, l, value)

im(1:l,:,:) = value;
im(:,1:l,:) = value;
im((end + 1) - (1:l), :,:) = value;
im(:, (end + 1) - (1:l),:) = value;