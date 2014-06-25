function launchMitoAnalysis(rootFolder, outputTemplate, summaryFile)

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
allIds = {};

for ii = 1:numel(dd)
    if ~(strcmp(dd(ii).name, '.') || strcmp(dd(ii).name, '..'))
        dirname = [rootFolder dd(ii).name];
        adir = dir([dirname '/*' extA]);
        mdir = dir([dirname '/*' extM]);
        
        if numel(adir) ~= numel(mdir)
            warning(['Found %d annotation files and %d mito files for' ...
                ' dir %s. Skipping.'], numel(adir), numel(mdir), dirname);
        else
            
            afiles = strcat([dirname '/'], sort({adir.name}));
            mfiles = strcat([dirname '/'], sort({mdir.name}));
            ids = repmat({dd(ii).name}, size(afiles));
            
            allAFiles = cat(2, allAFiles, afiles);
            allMFiles = cat(2, allMFiles, mfiles);
            allIds = cat(2, allIds, ids);
            
            % check for unused files.
            alldir = dir(dirname);
            fsel = [alldir.isdir] == 0;
            alldir = alldir(fsel);
            alldirfiles = {alldir.name};
            unusedfiles = setdiff(alldirfiles, cat(2, {adir.name}, {mdir.name}));
            if numel(unusedfiles) > 0
                fprintf('Did not process the following files for dir %s:\n',...
                    dirname);
                fprintf('%s\n', unusedfiles);
                fprintf('\n');
            end
        end
    end
end

terminalAndGliaMitoAnalysis(allAFiles, allMFiles, allIds, outputTemplate,...
    summaryFile);

end