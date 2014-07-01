function imstats_n = michelleAnalysis(annotationfiles_n, mitofiles_n,...
    animalid_n, outputTemplate, summaryfile)
% michelleAnalysis(annotationfiles, mitofiles,...
%    animalid, outputfile, summaryfile)
%
% Create terminal and glia mito analysis csv files
%   annotationfiles - a cell array containing the annotation image files
%   mitofiles       - a cell array containing the mitochondria image files
%   animalid        - a cell array of the animal ids
%   outputTemplate  - the template path to the output csv files
%                     should be in printf format with one string argument,
%                     like 'animalstuff_%s.csv'. Output files will be as:
%                     animalstuff_001.csv for the analysis file for the
%                                         animal with id 001, and
%                     animalstuff_001_obj.csv for the individual object
%                                         file for the same animal.
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
doMito_k = [true, true, false, false, false];

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

% Now, spoof the mito class
k = k + 1;
name_k{end + 1} = 'allMito';
imstats_n = spoofMito(imstats_n);

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
    sfname = sprintf(outputTemplate, animalid_a{i_a}); % stat file name
    ofname = sprintf(outputTemplate, ...
        sprintf('%s_obj', animalid_a{i_a})); % object file name
    
    m = numel(animal_stat_m);
    
    % Open the analysis file for writing and write a header to it
    fh = fopen(sfname, 'w');    
    if fh < 1
        error('Could not open %s for writing', sfname);
    end
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
    
    % Open the object file name for writing.
    oh = fopen(ofname, 'w');
    if oh < 1
        error('Could not open %s for writing', ofname);
    end
    % Write the object-specific header
    fprintf(oh, 'file, ');
    for i_k = 1:k
        fprintf(oh, '%s Area, %s Eccentricity, ',...
            name_k{i_k}, name_k{i_k});
    end
    fprintf(oh, '\n');
    
    % Now, write out the areas and eccenticities. Each object gets a
    % row in the column of its class. This code looks a little bit weird
    % because of that, just bare with us.
    
    
    for i_m = 1:m
        animal_stat = animal_stat_m(i_m);

        a = cell(1, k);
        ae = a;
        nk = zeros(1, k);
        
        for i_k = 1:k
            cname = name_k{i_k};
            a{i_k} = animal_stat.(cname).a;
            ae{i_k} = animal_stat.(cname).ae;
            nk(i_k) = numel(a{i_k});
        end
        
        p = max(nk);
        
        for i_p = 1:p
            fprintf(oh, '%s, ', animal_stat.afile);
            for i_k = 1:k                
                if nk(i_k) >= i_p
                    fprintf(oh, '%d, %f, ', a{i_k}(i_p), ae{i_k}(i_p));
                else
                    fprintf(oh, ', , ');
                end
            end
            fprintf(oh, '\n');
        end
    end
    
    fclose(oh);
    
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
function imstats_n = spoofMito(imstats_n)

global k name_k;

n = numel(imstats_n);

for i_n = 1:n
    a_mito = cell(1, k - 1);
    ae_mito = cell(1, k - 1);
    for i_k = 1:(k - 1)
        cname = name_k{i_k};
        a_mito{i_k} = imstats_n(i_n).(cname).m;
        ae_mito{i_k} = imstats_n(i_n).(cname).me;
    end
    imstats_n(i_n).allMito.a = [a_mito{:}];
    imstats_n(i_n).allMito.ae = [ae_mito{:}];
    imstats_n(i_n).allMito.m = [];
    imstats_n(i_n).allMito.me = [];    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imstats_n = collectStats(annotationfiles_n, mitofiles_n, ...
    animalid_n, n_files, k)

global minArea c_mito c_k name_k ids_n; 

% instantiate image struct array
% first, create a dummy struct
dummy = struct('afile', '', ... % annotation file
    'mfile', '', ...            % mitochondria file
    'animal', '', ...           % animal id (string)
    'animalidx', 0, ...         % animal id index (numeric)
    'totpix', 0);               % total pixels in the image
for i_k = 1:k
    cname = name_k{i_k};
    dummy.(cname).a = [];       % annotation object area
    dummy.(cname).ae = [];      % annotation object eccentricity
    dummy.(cname).m = [];       % object-associated mitochondria area
    dummy.(cname).me = [];      % --"-- eccentricity
end
    
imstats_n = repmat(dummy, size(annotationfiles_n));


parfor i_n = 1:n_files
    fprintf('Analyzing file %s\n', annotationfiles_n{i_n});
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
        
        [a, ae] = areaArray(a_mask);
        [m, me] = areaArray(m_mask);

        imstats_n(i_n).(cname).a = a;
        imstats_n(i_n).(cname).ae = ae;
        imstats_n(i_n).(cname).m = m;
        imstats_n(i_n).(cname).me = me;

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
function [a, ae] = areaArray(mask)
areaprop = regionprops(mask, 'Area', 'Eccentricity');
a = [areaprop.Area];
ae = [areaprop.Eccentricity];
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