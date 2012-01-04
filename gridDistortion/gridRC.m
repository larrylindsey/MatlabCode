function rc = gridRC(r, c)

[R C] = meshgrid(r, c);
rc = cat(2, R(:), C(:));
