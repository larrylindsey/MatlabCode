function michelleAnalysis(annotationfiles_n, mitofiles_n,...
    animalid_n, outputTemplate, summaryfile)
% michelleAnalysis(annotationfiles, mitofiles,...
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

fprintf('Begin analysis\n');

global minArea c_mito c_cap c_k name_k ids_n k strfields_f f doMito_k;


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
name_k = {'terminal', 'glia', 'capillary', 'background'};
strfields_f = {'n', 'fraction', 'avgsize'};
f = numel(strfields_f);

k = numel(name_k);

% If k_doMito[ii] is true, then do the mito analysis for the ii'th class
doMito_k = [true, true, false, false];

n_files = numel(annotationfiles_n);

% This line is confusing. Look at the output of help('unique') for guidance
[animalid_a, ~, ids_n] = unique(animalid_n);
a = numel(animalid_a);
% ids_a = 1:a

% Indices:
% n - file
% a - animal
% k - class
% f - struct fields

imstats_n = collectStats(annotationfiles_n, mitofiles_n, ...
    animalid_n, n_files, k);
% we should have
% all(ids_n == [imstats_n.animalidx])

% Create summary file, write the header
gh = fopen(summaryfile, 'w');
if gh < 1
    error('Could not open %s for writing', summaryfile);
end
fprintf(gh, 'Animal ID, ');
appendHeader(gh);

for i_a = 1:a
    animal_sel = ids_n == i_a;
    animal_stat_m = imstats_n(animal_sel);
    animal_summary_stat = computeSummary(animal_stat_m);
    fname = sprintf(outputTemplate, animalid_a{i_a});
    
    m = numel(animal_stat_m);
    
    fh = fopen(fname, 'w');
    
    if fh < 1
        error('Could not open %s for writing', fname);
    end
    
    % Write out the header
    fprintf(fh, 'Animal ID, Annotation File, Mito File, ');
    appendHeader(fh);
   
    % Now, write the data
    for i_m = 1:m
        animal_stat = animal_stat_m(i_m);
        
        fprintf(fh, '%s, %s, %s, ', animalid_a{i_a}, ...
            animal_stat.afile, animal_stat.mfile);
        
        appendData(fh, animal_stat);       
    end
    fclose(fh);
    
    fprintf(gh, '%s, ', animalid_a{i_a});
    appendData(gh, animal_summary_stat);
end

fclose(gh);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stat = computeSummary(stat_m)
global k name_k

stat.afile = {stat_m.afile};
stat.mfile = {stat_m.mfile};
stat.animal = stat_m(1).animal;
stat.totpix = sum([stat_m.totpix]);

for i_k = 1:k
    cname = name_k{i_k};
    c_cat = [stat_m.(cname)];
    stat.(cname).a = [c_cat.a];
    stat.(cname).m = [c_cat.m];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stat = computestat(sname, cname, str, type)
% sname - stat name
% cname - class name
% str - the struct
% type - 'a' or 'm' for annotation or mito

global minArea;

% population
pop = str.(cname).(type);
pop(pop < minArea) = [];

switch sname
    case 'n'
        stat = numel(pop);
    case 'fraction'
        % non capillary pixels
        noncappix = str.totpix - sum(str.capillary.a);
        stat = sum(pop) / noncappix;
    case 'avgsize'
        stat = mean(pop);
    otherwise
        error('don''t know how to compute stat %s', sname);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function appendData(hh, str)
global k name_k strfields_f f doMito_k;
for i_k = 1:k
    for i_f = 1:f
        fprintf(hh, '%g, ',...
            computestat(strfields_f{i_f}, name_k{i_k}, str, 'a'));
    end
    if doMito_k(i_k)
        for i_f = 1:f
            fprintf(hh, '%g, ',...
                computestat(strfields_f{i_f}, name_k{i_k}, str, 'm'));
        end
    end
end
fprintf(hh, '\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function appendHeader(hh)
global k f doMito_k name_k strfields_f;
for i_k = 1:k
    for i_f = 1:f
        fprintf(hh, '%s %s, ', name_k{i_k}, strfields_f{i_f});
    end
    if doMito_k(i_k)
        for i_f = 1:f
            fprintf(hh, '%s mito %s, ', name_k{i_k}, strfields_f{i_f});
        end
    end    
end
fprintf(hh, '\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imstats_n = collectStats(annotationfiles_n, mitofiles_n, ...
    animalid_n, n_files, k)

global minArea c_mito c_k name_k ids_n; 

% instantiate image struct array
imstats_n = repmat(struct, size(annotationfiles_n));


for i_n = 1:n_files
    % Annotation image
    im_a = im2single(imread(annotationfiles_n{i_n}));
    % Annotation image for mitochondria
    im_m = im2single(imread(mitofiles_n{i_n}));
    
    % Get the mitochondria mask from the mito image
    mito_mask = getMask(im_m, c_mito);
    % Get the capillary mask from the annotation image
    
    imstats_n(i_n).afile = annotationfiles_n{i_n};
    imstats_n(i_n).mfile = mitofiles_n{i_n};
    imstats_n(i_n).animal = animalid_n{i_n};
    imstats_n(i_n).animalidx = ids_n(i_n);
    
    imstats_n(i_n).totpix = size(im_a,1) * size(im_a,2);
    
    % Pixel analysis count
    npix = 0;
    
    % For each of the k classes
    for i_k = 1:k
        % Get the classes color
        c = c_k(i_k, :);
        cname = name_k{i_k};
        
        % Get the class-mask from the annotation image
        a_mask = getMask(im_a, c);
        m_mask = a_mask & mito_mask;
        
        imstats_n(i_n).(cname).a = areaArray(a_mask);
        imstats_n(i_n).(cname).m = areaArray(m_mask);
        
        npix = npix + sum(imstats_n(i_n).(cname).a);
    end
    
    if imstats_n(i_n).totpix - npix > minArea
        error(['For image file %s, missed some pixels. Analyzed %d, ', ...
            'expected %d'], annotationfiles_n{i_n}, npix, ...
            imstats_n(i_n).totpix);
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = areaArray(mask)
areaprop = regionprops(mask, 'Area');
a = [areaprop.Area];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = getMask(im, c)
% Creates a mask over im for a given color

% fprintf('Processing color %03d %03d %03d...', c(1) * 255, c(2) * 255, c(3)* 255);

% Start out everywhere true.
mask = true(size(im,1), size(im,2));
for ic = 1:numel(c)
    % Set the mask to false if the image's channel
    % doesn't match the corresponding component of
    % the color. This works for RGB or grayscale.
    mask = and(mask, abs(im(:,:,ic) - c(ic)) < (16 / 256));
end


end