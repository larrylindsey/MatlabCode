function [allAFiles, allMFiles] = getMichelleImages(rootFolder)

if nargin < 1
    rootFolder = './';
end

if rootFolder(end) ~= '/' && rootFolder(end) ~= '\'
    rootFolder = [rootFolder '/'];
end

% List the files in rootDir, then select only the ones that are directories
dd = dir(rootFolder);
dirsel = [dd.isdir] == 1;
dd = dd(dirsel);
extA = 'noMit.tif';
extM = 'wMit.tif';

if rootFolder(end) ~= '/' && rootFolder(end) ~= '\'
    rootFolder = [rootFolder '/'];
end

% List the files in rootDir, then select only the ones that are directories
dd = dir(rootFolder);
dirsel = [dd.isdir] == 1;
dd = dd(dirsel);

allAFiles = {};
allMFiles = {};

for ii = 1:numel(dd)
    if ~(strcmp(dd(ii).name, '.') || strcmp(dd(ii).name, '..'))
        dirname = [rootFolder dd(ii).name];
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