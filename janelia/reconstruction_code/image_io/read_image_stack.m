function [ imarray ] = read_image_stack( file_pattern, params )
% READ_IMAGE_STACK: read in a stack of images
% Multi-image TIFFs are supported, as are stacks of single images. The code
% makes the following assumptions:
%    1. If multiple images are to be loaded, they all have the same dimensions
%    and are all aligned.
%    2. If multiple files are to be read in, the filenames are in
%    lexicographical order, so that adjacent images in z correspond to adjacent
%    filenames as listed by the Matlab 'dir' function.
%    3. The images are either 8-bit or 16-bit B&W TIFFs
%    4. The file_pattern will match either a single (possibly multi-image)
%    TIFF file, or multiple single image files, but NOT multiple
%    multi-image files.
%
% INPUT:
%    - file_pattern: an expression interpretable by bash, possibly
%    including wild card characters, specifying the filename(s).
%    - params: OPTIONAL parameters for reading in the files, as follows:
%       * directory: the filesystem directory where the files are found.
%       * first_image: the first image to read, 1-indexed. 
%       * last_image: the last image to read, inclusive, 1-indexed.
%       * crop: a 2x2 matrix containing the first and last x coordinates on
%       the first row and the first and last y coordinates on the second
%       row (not supported for png)
%       * double: convert the intensity values to type double in [0, 1]
%       (flag)
%       * adjust: expand the pixel values in the image to take up the full
%       range of possible values (flag)
%       * histeq: perform uniform histogram equalization on the images
%       (flag)
%       * invert: take the negative (im_out = 1-im_in) of the input values
%       (flag)
%
% OUTPUT:
%    - imarray: a 3D array of stacked intensity values

if ~isfield(params, 'debug'),
    params.debug = 0;
end

% "intelligently" determine whether we are using a multi-image TIFF or
% multiple single-image TIFFs.
seps = strfind(file_pattern, '/');
if ~isempty(seps),
    directory_in_fn = [file_pattern(1:(seps(end)-1)) '/'];
    file_pattern = file_pattern((seps(end)+1):end);
else
    directory_in_fn = './';
end
if isfield(params, 'directory'),
    directory = [params.directory '/' directory_in_fn];
else
    directory = directory_in_fn;
end
filenames = dir([directory file_pattern]); % get the filenames
if params.debug,
    fprintf('The file pattern sought is: \n%s.\n', [directory file_pattern]);
    fprintf('The number of matching files is: %d.\n', numel(filenames));
end
image_info = imfinfo([directory filenames(1).name]);
if numel(filenames) == 1 && numel(image_info) >= 1,
    single_file_multi_image = true;
    num_images = numel(image_info);
elseif numel(filenames) >= 1 && numel(image_info) == 1,
    single_file_multi_image = false;
    num_images = numel(filenames);
else
    error(['Both the number of files and the number of images in the ' ...
        'first file are greater than 1.']);
end

if ~isfield(params, 'crop'),
    picXlen = image_info.Width;
    picYlen = image_info.Height;
    crop = { [1 picYlen], [1 picXlen] };
else
    crop = {params.crop(1,:), params.crop(2,:)};
    picXlen = crop{2}(2) - crop{2}(1) + 1;
    picYlen = crop{1}(2) - crop{1}(1) + 1;
end

%% Some sanity checks

if ~isfield(params, 'first_image'),
    params.('first_image') = 1;
end
if ~isfield(params, 'last_image'),
    params.('last_image') = num_images;
end
if params.first_image < 1 || params.first_image > num_images || ...
        params.last_image < 1 || params.last_image > num_images || ...
        params.last_image < params.first_image,
    error(['The first and last image specifications are invalid. ' ...
        'First image: %d; Second image: %d; Number of images: %d.'], ...
        params.first_image, params.last_image, num_images);
end
picZlen = params.last_image - params.first_image + 1;

%% Initialize imarray

switch image_info(1).BitDepth
    case 8
        imarray = zeros([picYlen picXlen picZlen], 'uint8');
        iscolor = false;
    case 16
        imarray = zeros([picYlen picXlen picZlen], 'uint16');
        iscolor = false;
    case 24
        imarray = zeros([picYlen picXlen picZlen 3], 'uint8');
        iscolor = true;
    otherwise
        error('Unknown image type.');
end

%% Read in the image data

if single_file_multi_image,
    for i=params.first_image:params.last_image,
        idx = i-params.first_image+1;
        if iscolor,
            imarray(:,:,idx,:) = ...
                imread([directory '/' file_pattern], i);
        else
            imarray(:,:,idx) = ...
                imread([directory '/' file_pattern], i, 'PixelRegion', crop);
        end
    end
else
    for i=params.first_image:params.last_image,
        idx = i-params.first_image+1;
        if iscolor,
            imarray(:,:,idx,:) = ...
                imread([directory '/' filenames(i).name]);
        else
            if strcmp(image_info(1).Format, 'png'),
                imarray(:,:,idx) = imread([directory '/' filenames(i).name]);

            else
                imarray(:,:,idx) = imread([directory '/' filenames(i).name], ...
                        'PixelRegion', crop);
            end
        end
    end
end

% Convert to double between 0 and 1
if isfield(params, 'double') && params.double,
    if iscolor,
        imarray = double(imarray)/(2^image_info(1).BitDepth/3-1);
    else
        imarray = double(imarray)/(2^image_info(1).BitDepth-1);
    end
end

% stretch to the full dynamic range of the data format
if isfield(params, 'adjust') && params.adjust,
    for iZ=1:size(imarray,3),
        if iscolor,
            pic = squeeze(imarray(:,:,iZ,:));
            imarray(:,:,iZ,:) = imadjust(pic,stretchlim(pic, 0));
        else
            pic = imarray(:,:,iZ);
            imarray(:,:,iZ) = imadjust(pic,stretchlim(pic, 0));
        end
    end
end

% equalize histogram to maximum entropy
if isfield(params, 'histeq') && params.histeq,
    for iZ=1:size(imarray,3),
        if iscolor,
            imarray(:,:,iZ,:) = histeq(squeeze(imarray(:,:,iZ,:)));
        else
            imarray(:,:,iZ) = histeq(imarray(:,:,iZ));
        end
    end
end

% take the negative of the image
if isfield(params, 'invert') && params.invert,
    imarray = max(imarray(:))-imarray;
end

end

