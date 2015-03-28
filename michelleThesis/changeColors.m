function changeColors(inDir, outDir)

oldColor = [255 255 0; 255  0    243] / 255;
newColor = .75 * [0 255 0; 255  0    243] / 255;

if nargin < 1
    inDir = './';
end

if inDir(end) ~= '/' && inDir(end) ~= '\'
    inDir = [inDir '/'];
end

% List the files in rootDir, then select only the ones that are directories
extA = 'noMit.tif';
extM = 'wMit.tif';

if inDir(end) ~= '/' && inDir(end) ~= '\'
    inDir = [inDir '/'];
end

% List the files in rootDir, then select only the ones that are directories
dd = dir(inDir);
dirsel = [dd.isdir] == 1;
dd = dd(dirsel);

allAFiles = {};
allMFiles = {};

for ii = 1:numel(dd)
    if ~(strcmp(dd(ii).name, '.') || strcmp(dd(ii).name, '..'))
        dirname = [inDir dd(ii).name];
        adir = dir([dirname '/*' extA]);
        mdir = dir([dirname '/*' extM]);
        
        if numel(adir) ~= numel(mdir)
            error(['Found %d annotation files and %d mito files for' ...
                ' dir %s.'], numel(adir), numel(mdir), dirname);
        else
            
            afiles = strcat([dirname '/'], sort({adir.name}));
            mfiles = strcat([dirname '/'], sort({mdir.name}));

            if ~isempty(afiles)               
                allAFiles = cat(2, allAFiles, afiles);
                allMFiles = cat(2, allMFiles, mfiles);
            end
        end
    end
end

allFiles = {allAFiles{:}, allMFiles{:}};

for i_f = 1:numel(allFiles)
    fqfn = allFiles{i_f};
    [folder, filename, ext] = fileparts(fqfn);
    if ~isempty(strfind(folder, inDir))
        outfolder = [outDir '/' folder((numel(inDir) + 1):end)];
    else
        outfolder = [outDir '/' folder];
    end
    system(['mkdir -p ' outfolder]);
    
    newFqfn = [outfolder '/' filename ext];
    
    fprintf('Changing %s to %s...\n', fqfn, newFqfn);
    
    im = imread(fqfn);
    im = replaceColor(im, oldColor, newColor);
    imwrite(im, newFqfn);
end

end
