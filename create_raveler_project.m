function control = create_raveler_project(control)
%

%folder, filter, outdir, control

%This code contains portions copied from a script written by Shiv
%Vitaladenuvi, at Janelia Farm, July 2010.

if nargin < 1
    control = defaultControl;
    return;
end

if ~isfield(control, 'output') || isempty(control.output)
    control.output = uigetdir(pwd);
    drawnow;
end

testrun = control.output(1) == 0;

%%% Begin Processing %%%

%last_boundary_energy = [];
links = [];
last_segment_map = [];
segment_offset = 0;
if isfield(control, 'debug')
    debug = control.debug;
else
    debug = false;
end

f = 1;
n = numel(control.files);

if testrun && (~isfield(control, 'index') || isempty(control.index))
    userinput = inputdlg('Enter Slice Index for Test', ...
        'Slice Index', 1, {'1'});
    control.index = str2double(userinput{1});    
end

if testrun
    f = control.index;
    n = control.index;
    debug = true;
    fprintf('Performing Test Run\n');
end

fprintf('Creating Raveler project for indices %d to %d\n', f, n);

control.section_ids = f:n;

for i_im = f:n
    fprintf('Iteration %d\n', i_im);
    control.index = i_im;
    % Get superpixel map, boundary energy
    fprintf('Boundary Energy\n');
    [boundary_energy control] = control.getBoundaryEnergy(control);
    fprintf('Super Pixels\n');
    [superpixel_map control] = control.getSuperPixels(control);

    fprintf('Segment map\n');
    % Create superpixel -> segment map
    [label_map, sp_to_seg] = ...
        compute_segmentation_from_superpixels_with_min_area_c(...
        superpixel_map, boundary_energy,...
        control.spseg.minimum_boundary_threshold,...
        control.spseg.minimum_area_threshold);%#ok

    % Calculate links
    sp_to_seg(:,2) = sp_to_seg(:,2) + segment_offset;
    segment_offset = max(sp_to_seg(:,2));

    superpixel_to_segment_maps{i_im} = sp_to_seg; %#ok

    segment_map = apply_mapping(superpixel_map, sp_to_seg);
    segment_map = double(remove_merged_boundaries_2D(uint32(segment_map)));
    segment_map(boundary_energy >...
        control.splink.boundary_energy_threshold) = 0;
    
    if ~isempty(last_segment_map)
        fprintf('Links\n');
        label_pairs_pixels = [segment_map(:), last_segment_map(:)];
        [label_pairs, overlap_area] = ...
            count_row_occurence(label_pairs_pixels);

        % include links
        links = [links; label_pairs, overlap_area]; %#ok
    end

    % Prepare next iteration
    % last_boundary_energy = boundary_energy;
    last_segment_map = segment_map;

    % Write image, superpixel map.
    if debug
        fprintf('View\n');
        close all;
        im = control.getData(control);
        plot_segment_boundaries(im, superpixel_map);
        drawnow;
    end
    
    if testrun
        return;        
    else
        fprintf('Write Images\n');
        write_image(control);
        write_spmap(superpixel_map, control);
    end
end

links = links(links(:,1)>0 & links(:,2)>0, :);

max_segment_id = segment_offset;

fprintf('building RAG ...\n');
rag = sparse([], [], [], max_segment_id, max_segment_id, 2*size(links, 1));
rag(sub2ind(size(rag), links(:,1), links(:,2))) = links(:,3);
rag(sub2ind(size(rag), links(:,2), links(:,1))) = links(:,3);

segment_to_body_map = get_connected_components(rag,...
    control.splink.area_overlap_threshold);

write_spsegmap(superpixel_to_segment_maps, control);

write_segbodymap(max_segment_id, segment_to_body_map, control);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [im control]= read_image(control)

if isfield(control, 'dataindex')
    lastindex = control.dataindex;
else
    lastindex = -1;
end

if lastindex == control.index
    im = control.dataImage;
else
    filestruct = control.files(control.index);

    imfile = filestruct.dfile;

    im = imread(imfile);
    if size(im, 3) > 1
        im = rgb2gray(im);
    end

    im = im2double(im);
    control.dataImage = im;
    control.dataindex = control.index;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_image(control)
im = control.getData(control);
check_for_dir(sprintf('%s/grayscale_maps/', control.output));
outfile = sprintf('%s/grayscale_maps/grayscale.%05d.png', ...
    control.output, control.index);

imwrite(im, outfile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_spmap(spmap, control)
check_for_dir(sprintf('%s/superpixel_maps/', control.output));
outfile = sprintf('%s/superpixel_maps/sp_map.%05d.png', control.output,...
    control.index);
imwrite(spmap, outfile, 'BitDepth', 16);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_spsegmap(superpixel_to_segment_maps, control)
fprintf('writing superpixel to segment map ...\n');
fout = fopen([control.output, '/superpixel_to_segment_map.txt'], 'wt');
for z = control.section_ids
    fprintf(fout, '%d\t%d\t%d\n', ...
        [repmat(z, [1 1+ size(superpixel_to_segment_maps{z}, 1)]); ...
        zeros(2, 1), superpixel_to_segment_maps{z}']);
end
fclose(fout);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_segbodymap(max_segment_id, segment_to_body_map, control)
fprintf('writing segment to body map ...\n');
fout = fopen([control.output, '/segment_to_body_map.txt'], 'wt');
fprintf(fout, '%d\t%d\n', [zeros(2,1),...
    [1:max_segment_id; segment_to_body_map']]);
fclose(fout);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [be control] = readPredictionFile(control)

if isfield(control, 'beindex')
    lastindex = control.beindex;
else
    lastindex = -1;
end

if control.index == lastindex
    be = control.boundaryEnergyImage;
else
    filestruct = control.files(control.index);
    befile = filestruct.pfile;

    be = imread(befile);

    be = im2double(be);
    
    control.boundaryEnergyImage = be;
    
    control.beindex = control.index;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [segmap control] = defaultSegmentation(control)

[im control] = control.getData(control);
[be control] = control.getBoundaryEnergy(control);

segmap = better_segment(be, im);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [segmap control] = readSegmentation(control)
index = control.index;
segmap = im2double(imread(control.files(index).sfile));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function control = defaultControl
control.input = '';
control.output = '';
control.filter = '';
control.getBoundaryEnergy = @readPredictionFile;
control.getSuperPixels = @defaultSegmentation;
control.getData = @read_image;
control.spseg.minimum_boundary_threshold = 0.1;
control.spseg.minimum_area_threshold = 2500;
control.splink.boundary_energy_threshold = 0.3;
control.splink.area_overlap_threshold = 100;
control.files.dfile = '';
control.files.pfile = '';
control.altsp = @readSegmentation;
end