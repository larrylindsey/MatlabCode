function [stat_str, masks] = imagesegmentstats(imagefile, outputfile)
% stat_str imageSegmentStats(imagefile, outputfile)
% Calculates segment stats on a given image file
%    imagefile - the path to the image file to analyze
%                pass as empty to run the file dialog
%    outputfile - the path to the output csv file
%                 if empty, no file is generated
%    stat_str - a struct array containing the blob stats

% global COLOR_RGB_DIST;
% COLOR_RGB_DIST = 3;

% Predefined colors
im_kc = [255     0     0;
    184    0     28;
    26    93     0;
    43    93     0;
    255   255     0;
    0     0   255;
    146     6   181;
    96     31   151;
    255     0   243;
    0   255   255;
    255   153     0;
    0     0     0] / 255;

if isempty(imagefile)
    [fname, pname] = uigetfile();
    
    if pname(end) ~= '/'
        pname = [pname '/'];
    end
    imagefile = [pname fname];
end

% Read the image file
if ischar(imagefile)
    fprintf('Reading %s...', imagefile);
    im = imread(imagefile);
    imsz = size(im);
    totPix = imsz(1) * imsz(2);
    fprintf('Done. Read %d pixels\n', totPix);
else
    im = imagefile;
    imsz = size(im);
    totPix = imsz(1) * imsz(2);
end

im = im2single(im);

% % Get a list of the unique colors in the image
% im_kc = reshape(im, [prod(imsz(1:2)) imsz(3)]);
% im_kc = unique(im_kc, 'rows');

% k is the number of unique colors
k = size(im_kc, 1);

%stat_str = {struct};
%stat_str = repmat(stat_str, [1 k]);
stat_str = repmat(struct,0);

masks = false(size(im,1), size(im,2), k);

% iterate through the colors
for ik = 1:k
    fprintf('Color %d of %d\n', ik, k);
    % c is a 1 or 3 element vector containing the current color
    c = im_kc(ik,:);
    
    fprintf('Generating mask...\n');
    % mask is a logical image that is true where im is the current
    % color, and false where it is not.
    mask = getMask(im, c);
    %   extracted_mask = or(extracted_mask, mask);
    masks(:,:,ik) = mask;
    
    fprintf('Collecting region properties...\n');
    % str_ik is a struct array with the field Area.
    % There is one element in the array per connected-component
    % in the image, and Area corresponds to its area, in pixels
    str_ik = regionprops(mask, 'Area');
    
    fprintf('Collating statistics...\n');
    % Black or white is expected to be background
    % In this case, clump the areas into a single sum
    if all(c == 0) || all(c == max(im_kc(:)))
        sumArea = sum([str_ik.Area]);
        str_ik = str_ik(1);
        str_ik.Area = sumArea;
    end
    
    % Add fields for group number, color and fraction to each struct
    for is = 1:numel(str_ik)
        str_ik(is).Color = c;
        str_ik(is).idx = ik;
        str_ik(is).Fraction = str_ik(is).Area / totPix;
    end
    
    % Concatenate the result into the return value
    if (numel(str_ik) > 0)
        stat_str = cat(1, stat_str, str_ik);
    end
    
    fprintf('Done.\n');
end

% stat_str = cat(1, stat_str{:});

% If we're given an output file
if nargin > 1 && ~isempty(outputfile)
    % open the file
    f = fopen(outputfile, 'w');
    % f is the file handle. f < 1 indicates an error of some kind
    if f < 1
        fprintf('Could not open file %s for writing\n', outputfile);
    else
        % Write the header row
        fprintf(f, 'Group ID, Color, Pixel Count, Fraction of Image\n');
        % Then write the data rows
        for is = 1:numel(stat_str)
            str_i = stat_str(is);
            fprintf(f, '%d, %s, %d, %g\n', str_i.idx, hexColor(str_i.Color * 255), ...
                str_i.Area, str_i.Fraction);
        end
        % close and flush the file. We're done.
        fclose(f);
    end
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = getMask(im, c)
% Creates a mask over im for a given color

fprintf('Processing color %03d %03d %03d...', c(1) * 255, c(2) * 255, c(3)* 255);

% Start out everywhere true.
mask = true(size(im,1), size(im,2));
for ic = 1:numel(c)
    % Set the mask to false if the image's channel
    % doesn't match the corresponding component of
    % the color. This works for RGB or grayscale.
    mask = and(mask, abs(im(:,:,ic) - c(ic)) < (16 / 256));
end

fprintf('done.\n');

%subplot(1, 2, 1); imshow(imresize(im, 1/24));
%subplot(1, 2, 2); imshow(imresize(mask, 1/24));
%drawnow;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hex = hexColor(c)
% Create a color string that can be used in Excel

% If we get only one component, set all three to
% be the same. We use < 3 rather than == 1 for
% robustness. Plus, it looks like a heart. Kind of.
if numel(c) < 3
    repmat(c, [1 3]);
end



hex = sprintf('#%02X%02X%02X', c(1), c(2), c(3));

end