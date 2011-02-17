function merge_images(leftdir, rightdir, mergedir, fmt)

if nargin < 4
    if nargin < 3
        if nargin < 2
            if nargin < 1
                leftdir = uigetdir(pwd, 'Left Panel Images');
                if leftdir(1) == 0
                    return
                end
            end
            rightdir = uigetdir(leftdir, 'Right Panel Images');
            if rightdir(1) == 0
                return
            end
        end
        mergedir = uigetdir(rightdir, 'Output Directory');
        if mergedir(1) == 0
            return
        end
    end
    fmt = questdlg('Output which image format?',...
        'Output Image Format', 'jpg','tiff','jpg');
else
    if not(or(strncmpi(fmt,'jpg',3),strncmpi(fmt,'tif',3)))
        error(['%s is not a valid image format.  Valid image formats' ...
            ' are jpg or tiff'], fmt);
    end
end

rlist = dir(rightdir);
llist = dir(leftdir);

rlist = rlist(3:end);
llist = llist(3:end);

nfiles = min(length(llist), length(rlist));

if length(rlist) ~= length(llist)
    res = questdlg('Directories contain different numbers of files.', ...
        'Uh Oh', 'Continue', 'Cancel', 'Cancel');
    if strcmp(res, 'Cancel');
        return;
    end
end

logfid = fopen([mergedir '/log.csv'], 'w');
fprintf(logfid, ['Left Image,Right Image,Merged Image,Confidence,',...
    'X Shift,Y Shift\n']);

wh = waitbar(0, 'Initializing...');
whstr = handle2struct(wh);
pp = get(wh, 'Position');
pp(end) = 100;
set(wh, 'Position', pp);

set(whstr.children.children(1).handle,'Interpreter','none');

for ilist = 1:nfiles
    lname = llist(ilist).name;
    rname = rlist(ilist).name;
    mergename = [lname ' merged ' rname];
    mergename(logical(mergename=='.')) = '_';
    mergename = [mergename '.' fmt];

    
    waitbar(ilist / nfiles, wh, ...
        sprintf('Merging %s and %s into\n %s.\nClose this window to cancel',...
        lname, rname, mergename));
    drawnow;
    
    imlt = imread([leftdir '/' lname]);
    imlt = mean(imlt, 3);
    imrt = imread([rightdir '/' rname]);
    imrt = mean(imrt, 3);
    [immerge conf xs ys] = montage_images_displaced(imlt, imrt);
    imwrite(immerge / max(immerge(:)), [mergedir '/' mergename]);
    clear immerge;
    logstr = sprintf('%s,%s,%s,%g,%g,%g\n',...
        lname, rname, mergename, conf, xs, ys);
    fprintf(logfid, logstr);
    
    if ~ishandle(wh)
        return;
    end
end
fclose(logfid);
end