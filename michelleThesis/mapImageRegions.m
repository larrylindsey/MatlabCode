function files = mapImageRegions(rootFolder)

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

files = {};

for ii = 1:numel(dd)
    files{ii} = {};
    if ~(strcmp(dd(ii).name, '.') || strcmp(dd(ii).name, '..'))
        dirname = [rootFolder dd(ii).name];
        tifdir = dir([dirname '/*.tif']);
        tifs = {tifdir.name};
        files = cat(2, files{:}, strcat([dirname '/'], tifs));
    end
end

end