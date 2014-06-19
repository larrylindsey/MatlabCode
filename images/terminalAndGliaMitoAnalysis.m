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


minArea = 256; % corresponds to blob roughly 16 x 16 pix.

% Predefined colors

% Mito color (blue)
c_mito = [0 0 1];

% Colors among k number of classes
c_k = [255  0	0;
       255	255	0;
       255  0   243;
       0    0   0] / 255;

% class names
kname = {'terminal', 'glia', 'capillary', 'background'};
strfields = {'n', 'fraction', 'avgsize'};

k = numel(kname);

% If k_doMito[ii] is true, then do the mito analysis for the ii'th class
k_doMito = [true, true, false, false];

n_files = numel(annotationfiles);

stat_str = repmat(areaStatStr, [n_files, k]);
mito_str = stat_str;

[uanimalid, ~, ids] = unique(animalid);

for ii = 1:n_files
    im_a = im2single(imread(annotationfiles{ii}));
    im_m = im2single(imread(mitofiles{ii}));
    
    mito_mask = getMask(im_m, c_mito);
    
    id = ids(ii);
    
    for ik = 1:k
        c = c_k(ik, :);
        mask = getMask(im_a, c);
        
        stat_str(ii, ik) = areaStatStr(mask, minArea, id, ...
            annotationfiles{ii});        
        
        if k_doMito(ik)            
            mito_str(ii, ik) = areaStatStr(and(mask, mito_mask), minArea,...
                id, mitofiles{ii});
        end
        
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
    for ik = 1:k
        for ff = 1:numel(strfields)
            fprintf(g, '%g, ', ...
                mean([animal_stat_str(:, ik).(strfields{ff})]));
        end
        if k_doMito(ik)
            for ff = 1:numel(strfields)
                fprintf(g, '%g, ', ...
                    mean([animal_stat_str(:, ik).(strfields{ff})]));
            end
        end        
    end
    fprintf(g, '\n');
end

fclose(g);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = areaStatStr(mask, minArea, id, fname)

if nargin < 1 || isempty(mask)
    str = struct('n', 0, 'fraction', 0, 'avgsize', 0, 'id', -1, ...
        'filename', '');
    return
end

totPix = size(mask, 1) * size(mask, 2);
rpstr = regionprops(mask, 'Area');
area = [rpstr.Area];
sel = area < minArea;
area(sel) = [];

str.n = numel(area);
str.fraction = sum(area) / totPix;
str.avgsize = mean(area);
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
