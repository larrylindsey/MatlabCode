function mergeImages(files, prefix, type)

if nargin < 3
    type = 'png';
    if nargin < 2
        prefix = 'merged_image_';
    end
end

im{1} = adjustContrast(imread(files{1}));

n = numel(files);
w = floor(log(n) / log(10)) + 1;
subsStr = sprintf('%%0%gg', w);

for i_f = 2:numel(files)
    im{i_f} = adjustContrast(imread(files{i_f}));
    imCurr = im{i_f};
    imLast = im{i_f - 1};
    imMerge = (imCurr + imLast) / 2;
    
    fileName = [prefix sprintf(subsStr, i_f) '.' type];
    imwrite(imMerge, fileName);
    
end

end