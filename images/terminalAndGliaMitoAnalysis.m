    function terminalAndGliaMitoAnalysis(annotationfiles, mitofiles,...
    animalid, outputTemplate, summaryfile)
% terminalAndGliaMitoAnalysis(annotationfiles, mitofiles,...
%    animalid, outputfile, summaryfile)
%
% Create terminal and glia mito analysis csv files
%   annotationfiles - a cell array containing the annotation image files
%   mitofiles       - a cell array containing the mitochondria image files
%   animalid        - a cell array of the animal ids
%   outputTemplate  - the template path to the output csv file
%                     should be in printf format with one string argument,
%                     like 'animalstuff_%s.csv'
%   summaryfile     - the path the summary csv file
%
% annotationfiles, mitofiles, and animalid must all have the same number of
% elements.

%%%% Setup Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minArea = 256; % corresponds to blob roughly 16 x 16 pix.

% Predefined colors

% Mito color (blue)
c_mito = [0 0 1];

% Capillary color (purple / magenta)
c_cap = [255  0   243] / 255;

% Colors among k number of classes
c_k = [255  0	0;
       255	255	0;       
       255  0    243;
       0    0   0] / 255;

% class names
kname = {'terminal', 'glia', 'capillary', 'background'};
strfields = {'n', 'fraction', 'avgsize'};

k = numel(kname);

% If k_doMito[ii] is true, then do the mito analysis for the ii'th class
k_doMito = [true, true, false, false];

n_files = numel(annotationfiles);

% Instantiate the stat_str struct
stat_str = repmat(areaStatStr, [n_files, k]);
% and the mito_str struct
mito_str = stat_str;

% This line is confusing. Look at the output of help('unique') for guidance
[uanimalid, ~, ids] = unique(animalid);

%%%% Collate Data Structures from Image Files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:n_files
    % Annotation image
    im_a = im2single(imread(annotationfiles{ii}));
    % Annotation image for mitochondria
    im_m = im2single(imread(mitofiles{ii}));
    
    % Get the mitochondria mask from the mito image
    mito_mask = getMask(im_m, c_mito);    
    % Get the capillary mask from the annotation image
    cap_mask = getMask(im_a, c_cap);
    % Now we know how many pixels the capillary takes up
    capPix = sum(cap_mask(:));
    
    if capPix == 0
        warning('Found zero capillary pixels');
    end
    
    totPix = size(im_a, 1) * size(im_a, 2) - capPix;
    npix = zeros(k, 1);
    
    % Store the current animal id index
    id = ids(ii);
    
    % For each of the k classes
    for ik = 1:k
        % Get the classes color
        c = c_k(ik, :);
        % Get the class-mask from the annotation image
        mask = getMask(im_a, c);
        
        % Compute the statistic struct for this class
        [stat_str(ii, ik), npix(ik)] = areaStatStr(mask, minArea, id, ...
            annotationfiles{ii}, totPix);        
        
        % If we want the corresponding mitochondria stats, get those here
        % We compute over the pixel-wise and of the class mask and the mito
        % mask.
        if k_doMito(ik)            
            mito_str(ii, ik) = areaStatStr(and(mask, mito_mask), minArea,...
                id, mitofiles{ii}, totPix);
        end
        
    end
    
    if sum(npix) ~= size(im_a, 1) * size(im_a, 2)
        error(['For image file %s, missed some pixels. Analyzed %d, ', ...
            'expected %d'], annotationfiles{ii}, sum(npix), ...
            size(im_a, 1) * size(im_a, 2));
    end
    
end

strid = [stat_str(:,1).id];

g = fopen(summaryfile, 'w');

if g < 1
    error('Could not open %s for writing', filename);
end

fprintf(g, 'Animal ID, ');
for ik = 1:k
    for ff = 1:numel(strfields)
        fprintf(g, '%s %s, ', kname{ik}, strfields{ff});
    end
    if k_doMito(ik)
        for ff = 1:numel(strfields)
            fprintf(g, '%s mito %s, ', kname{ik}, strfields{ff});
        end
    end    
end
fprintf(g, '\n');

% Write out animal files
for ii = 1:numel(ids)
    id = ids(ii);
    sel = find(strid == id);
    animal_stat_str = stat_str(sel,:);
    animal_mito_str = mito_str(sel,:);
    
    n = numel(sel);
    
    fname = sprintf(outputTemplate, uanimalid{id});
    
    f = fopen(fname, 'w');
    
    if f < 1
        error('Could not open %s for writing', fname);
    else
        % Write out the header
        fprintf(f, 'Animal ID, Annotation File, Mito File, ');
        for ik = 1:k
            for ff = 1:numel(strfields)
                fprintf(f, '%s %s, ', kname{ik}, strfields{ff});
            end
            if k_doMito(ik)
                for ff = 1:numel(strfields)
                    fprintf(f, '%s mito %s, ', kname{ik}, strfields{ff});
                end
            end            
        end
        fprintf(f, '\n');
        
        for jj = 1:n
            fprintf(f, '%s, %s, %s, ', uanimalid{id},...
                animal_stat_str(jj, 1).filename, ...
                animal_mito_str(jj, 1).filename);
            for ik = 1:k
                for ff = 1:numel(strfields)
                    fprintf(f, '%g, ', ...
                        animal_stat_str(jj, ik).(strfields{ff}));
                end
                if k_doMito(ik)
                    for ff = 1:numel(strfields)
                        fprintf(f, '%g, ', ...
                            animal_mito_str(jj, ik).(strfields{ff}));
                    end
                end                
            end
            fprintf(f, '\n');
        end
        
    end
    
    fclose(f);
    
    fprintf(g, '%s, ', uanimalid{id});
    % ik is class index, k number of classes
    for ik = 1:k
        for ff = 1:numel(strfields)
            aggregate_stat = ...
                aggregate(animal_stat_str(:,ik), strfields{ff});
            fprintf(g, '%g, ', aggregate_stat);
        end
        if k_doMito(ik)
            for ff = 1:numel(strfields)
                aggregate_stat = ...
                    aggregate(animal_mito_str(:,ik), strfields{ff});
                fprintf(g, '%g, ', aggregate_stat);
            end
        end        
    end
    fprintf(g, '\n');
end

fclose(g);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = aggregate(str, fname)
% str - an animal stat struct array, one struct per image
% fname - the field name to aggregate. each field has its own rule
switch fname
    case 'n'
        % n aggregated by summation over all images
        a = sum([str(:, ik).n]);
    case 'fraction'
        % fraction aggregated by summing npix across all images, then 
        % dividing by the summation of totpix.
        % this gives the total fraction over all analyzed non-capillary
        % pixels.
        a = sum([str(:, ik).fraction] .* [str(:, ik).totpix]) /...
            sum([str(:, ik).totpix]);
    case 'avgsize'
        % avgsize aggregated by summing totsize, found by image-wise n *
        % avgsize, then dividing the sum by the n across all images.
        a = sum([str(:, ik).avgsize] .* [str(:, ik).n]) / ...
            sum([str(:, ik).n]);
    otherwise
        % As of this writing, we know of only three output struct fields.
        % If we encounter an unexpected one, aggregate to a not-a-number
        % (NaN) value.
        a = nan;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [str, npix] = areaStatStr(mask, minArea, id, fname, totPix)

if nargin < 1 || isempty(mask)
    str = struct('n', 0, 'fraction', 0, 'avgsize', 0, 'totpix', 0, ...
        'id', -1, 'filename', '');
    return
end

rpstr = regionprops(mask, 'Area');
area = [rpstr.Area];
npix = sum(area);
sel = area < minArea;
area(sel) = [];

str.n = numel(area);
str.fraction = sum(area) / totPix;
str.avgsize = mean(area);
str.totpix = totPix;
str.id = id;
str.filename = fname;

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
