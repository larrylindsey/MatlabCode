function write_image_stack(stack, params)
% WRITE_IMAGE_STACK: Write a stack of images to disk, either as a series of
%    8-bit TIFFs or PNGs, or as a single, multi-image TIFF. The series can
%    be either grayscale (size [m n p]) or RGB color (size [m n 3 p])
% INPUT:
%    - stack: the stack, a 3-dimensional matrix
%    - params: the writing parameters:
%        * filename: the filename for the output (setting this field
%        implies single_file, and is incompatible with prefix, padding, suffix)
%        * prefix, padding, suffix: the parameters that determine the output
%        filename(s). As an example, if prefix == 'image', padding == 4, and 
%        suffix == 'out', then if we are writing a stack of files, from 20 to 
%        49, then the filename will be image.0020.0049.out.tif, or if we are 
%        writing multiple files, then file 25 will be called 
%        image.0025.out.tif. Use padding == 1 if you don't want any padding.
%        * directory: the file system directory in which to write the files
%        * type: the type of image file, e.g. 'tif' or 'png'. Note that
%        single_file only works for 'tif'.
%        * first_filenum: the first number in the series. (stack(:,:,1) will be
%        written as the first_filenum^th file)
%        * first_slice: the first slice in the image stack to be written to
%        disk
%        * num_slices: the last number in the series.
%        (stack(:,:,first_slice+num_slices-1) will be written as the last file
%        * single_file: Whether the stack should be written as multiple, 
%        sequential TIFFs or a single, multi-image TIFF
%        * compression: should the files be compressed (directly passed to 
%        Matlab's imwrite function)
%        * invert: should the image files be inverted? (img_out=1-img_in)
%        * adjust: enforce that the image data is in [0,1] (default is
%        true)
%        * histeq: perform histogram equalization before writing

    %% Preprocessing of the parameters and default value assignments
    if isfield(params, 'invert') && params.invert,
        stack = 1-stack;
    end
    
    if isfield(params, 'compression'),
        compression = params.compression;
    else
        compression = 'none';
    end

    if ~isfield(params, 'single_file'),
        params.single_file = false;
    end

    if isfield(params, 'filename'),
        params.single_file = true;
        seps = strfind(params.filename, '/');
        if ~isempty(seps),
            directory_in_fn = [params.filename(1:(seps(end)-1)) '/'];
            params.filename = params.filename((seps(end)+1):end);
        else
            directory_in_fn = '';
        end
    else
        directory_in_fn = '';
    end
    
    if isfield(params, 'filename') && isfield(params, 'prefix'),
        error('Only filename or (prefix,padding,suffix) can be specified.');
    end
    
    if isfield(params, 'padding'),
        format_str = sprintf('%%0%dd', params.padding);
    end
    if ~isfield(params, 'type'),
        params.type = 'tif';
    end
    if ~isfield(params, 'adjust'),
        params.adjust = true;
    end
    if ~isfield(params, 'histeq'),
        params.histeq = false;
    end
    
    if isfield(params, 'directory'),
        directory = [params.directory '/' directory_in_fn];
    else
        directory = directory_in_fn;
    end
    
    if ~isfield(params, 'first_slice'),
        params.first_slice = 1;
    end
    
    if ~isfield(params, 'num_slices'),
        params.num_slices = size(stack, ndims(stack));
    end
    
    if isfield(params, 'first_filenum'),
        first_filenum = params.first_filenum;
        last_filenum = first_filenum+params.num_slices-1;
    else
        first_filenum = 0;
        last_filenum = first_filenum+params.num_slices-1;
    end
    
    if isfield(params, 'prefix') && isfield(params, 'suffix'),
        prefix = params.prefix;
        suffix = params.suffix;
    end
    
    s = size(stack);
 
    % if we are using 'double' for the image data, rescale if necessary to
    % [0, 1] before converting to uint8 for output.    
    if isa(stack, 'double'),
        if min(stack(:)) < 0,
            stack = stack - min(stack(:));
        end
        if max(stack(:)) > 1,
            stack = stack / max(stack(:));
        end
    end
    
    if ndims(stack) == 3,
        color_image = false;
        stack = reshape(stack, [s(1) s(2) 1 s(3)]);
    elseif ndims(stack) == 4,
        color_image = true;
    end
    
    s = size(stack);
    if params.adjust,
        stackp = permute(stack, [1 2 4 3]);
        stackp = reshape(stackp, [s(1) s(2)*s(4) s(3)]);
        stackp = reshape(imadjust(stackp, stretchlim(stackp, 0)), ...
                                 [s(1) s(2) s(4) s(3)]);
        stack = permute(stackp, [1 2 4 3]);
        clear stackp;
    end
    
    if params.histeq,
        stackp = permute(stack, [1 2 4 3]);
        stackp = reshape(stackp, [s(1) s(2)*s(4) s(3)]);
        stackp = reshape(histeq(stackp), [s(1) s(2) s(4) s(3)]);
        stack = permute(stackp, [1 2 4 3]);
        clear stackp;
    end
    
    if isa(stack, 'double'),
        stack = uint8(255*stack);
    end


    %% Writing of the files
    if params.single_file,
        if isfield(params, 'filename'),
            filename = params.filename;
        else
            filename = sprintf(...
                [prefix '.' format_str '.' format_str '.' suffix '.' ...
                params.type], first_filenum, last_filenum);
        end
        if ~color_image, % B&W image
            cur_image = squeeze(stack(:,:,params.first_slice));
        else % color image
            cur_image = squeeze(stack(:,:,:,params.first_slice));
        end
        imwrite(cur_image, [directory filename], params.type, 'Compression', ...
            compression);
        for filenum=(first_filenum+1):last_filenum,
            if ~color_image, % B&W image
                cur_image = squeeze( ...
                        stack(:,:,filenum-first_filenum+params.first_slice));
            else % color image
                cur_image = squeeze( ...
                        stack(:,:,:,filenum-first_filenum+params.first_slice));
            end
            imwrite(cur_image, [directory filename], params.type, ...
                'WriteMode', 'append', 'Compression', compression);
        end
    else
        for filenum=first_filenum:last_filenum,
            filename = [prefix '.' sprintf(format_str, filenum) ...
                                                    suffix '.' params.type];
            if color_image,
                cur_image = squeeze( ...
                        stack(:,:,:,filenum-first_filenum+params.first_slice));
            else
                cur_image = squeeze( ...
                        stack(:,:,filenum-first_filenum+params.first_slice));
            end
            imwrite(cur_image, [directory filename], params.type, ...
                'Compression', compression);
        end
    end
end
