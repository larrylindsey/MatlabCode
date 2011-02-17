function im = unsplitImage(outputImage)

[path prename ext] = fileparts(outputImage);

if strcmp(path, '')
    path = pwd;
end

% Calculate number of rows and columns
r = 1;
c = 1;
while exist(sprintf('%s/r%d_c%d', path, r, c), 'dir') > 0
    r = r + 1;
end
r = r - 1;

while exist(sprintf('%s/r%d_c%d', path, r, c), 'dir') > 0
    c = c + 1;
end
c = c - 1;

rowchunk = cell(1, r);
for i_r = 1:r
    colchunk = cell(1, c);
    for i_c = 1:c
        fname = sprintf('%s/r%d_c%d/%s_row%d_col%d%s', ...
            path, i_r, i_c, prename, i_r, i_c, ...
            ext);
        colchunk{i_c} = imread(fname);
        
    end
    rowchunk{i_r} = cat(2, colchunk{:});
end


im = cat(1, rowchunk{:});

if nargout < 1
    imwrite(im, outputImage);
end