function writeReconstruct(serdoc, secdoc, name)

if nargin < 3
    name = secdoc(1).name;
    idot = find(name == '.', 1, 'last');
    name = name(1:(idot-1));
end

serfid = fopen(sprintf('%s.ser', name), 'r');
if serfid > 0
    answer = '';
    while ~strcmpi(answer, 'yes') 
        answer = input(sprintf(...
            'Reconstruct series %s exists. Overwrite? (yes/no) ',...
            name), 's');
        if strcmpi(answer, 'no') || strcmpi(answer, 'n')
            return;
        end
    end
    fclose(serfid);
end

serfid = fopen(sprintf('%s.ser', name), 'w');
writeSer(serdoc, serfid);
fclose(serfid);

for index = 1:numel(secdoc)
    secfid = fopen(sprintf('%s.%d', name, secdoc(index).index), 'w');
    writeSec(secdoc(index), secfid);
    fclose(secfid);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = fieldToText(fname, s)

content = s.(fname);

if ischar(content)
    contentText = content;
elseif islogical(content)
    if content
        contentText = 'true';
    else
        contentText = 'false';
    end
elseif min(size(content)) == 1
    contentText = '';
    
    for ic = 1:numel(content)
        if abs(content(ic) - floor(content(ic))) == 0
            cstr = '%s %d';
        else
            cstr = '%s %g';
        end
        contentText = sprintf(cstr, contentText, content(ic));
    end
    contentText = stripOutterSpaces(contentText);
else
    contentText = '';
    for r = 1:size(content, 1)
        for c = 1:size(content, 2)
            contentText = sprintf('%s %g', contentText, content(r, c));
        end
        contentText = sprintf('%s,\n\t', contentText);
    end
end

text = sprintf('%s="%s"', fname, contentText);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeSer(serdoc, fid)

fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<!DOCTYPE Series SYSTEM "series.dtd">\n\n');
fprintf(fid, '<Series ');
names = fieldnames(serdoc);
for n = 1:numel(names)
    currfn = names{n};
    if strcmpi(currfn, 'Contour')
        fprintf(fid, '>\n');
        contours = serdoc.Contour;
        for ic = 1:numel(contours)
            currContour = contours(ic);
            fprintf(fid, '<Contour ');
            cnames = fieldnames(currContour);
            for icn = 1:numel(cnames)
                currcfn = cnames{icn};
                fprintf(fid, '%s\n', fieldToText(currcfn, ...
                    currContour));
            end
            fprintf(fid, '/>\n');
        end
    elseif strcmpi(currfn, 'ZContour')
        contours = serdoc.ZContour;
        for ic = 1:numel(contours)
            currContour = contours(ic);
            fprintf(fid, '<ZContour ');
            cnames = fieldnames(currContour);
            for icn = 1:numel(cnames)
                currcfn = cnames{icn};
                fprintf(fid, '%s\n', fieldToText(currcfn, ...
                    currContour));
            end
            fprintf(fid, '/>\n');
        end
    elseif ~strcmpi(currfn, 'pre')
        fprintf(fid, '%s\n', fieldToText(currfn, serdoc));
    end
end
fprintf(fid, '</Series>\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeSec(secdoc, fid)
section = secdoc.section;
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<!DOCTYPE Series SYSTEM "series.dtd">\n\n');
fprintf(fid, '<Section %s %s %s >\n', ...
    fieldToText('index', secdoc), fieldToText('thickness', section), ...
    fieldToText('alignLocked', section));

transforms = section.Transform;

for it = 1:numel(transforms)
    trans = transforms(it);
    fprintf(fid, '<Transform %s\n%s\n%s>\n', ...
        fieldToText('dim', trans), fieldToText('xcoef', trans), ...
        fieldToText('ycoef', trans));
    images = trans.Image;
    contours = trans.Contour;
    
    for ii = 1:numel(images)
        im = images(ii);
        fprintf(fid, '<Image %s %s %s %s %s %s %s />\n', ...
            fieldToText('mag', im), ...
            fieldToText('contrast', im),...
            fieldToText('brightness', im),...
            fieldToText('red', im),...
            fieldToText('green', im),...
            fieldToText('blue', im),...
            fieldToText('src', im));
    end
    
    for ic = 1:numel(contours)
        con = contours(ic);
        fprintf(fid, '<Contour %s %s %s %s %s %s %s %s />\n', ...
            fieldToText('name', con),...
            fieldToText('hidden', con),...
            fieldToText('closed', con),...
            fieldToText('simplified', con),...
            fieldToText('border', con),...
            fieldToText('fill', con),...
            fieldToText('mode', con),...
            fieldToText('points', con));
    end
    
    fprintf(fid, '</Transform>\n');
end

fprintf(fid, '</Section>\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = stripOutterSpaces(text)

fi = find(text ~= ' ', 1, 'first');
li = find(text ~= ' ', 1, 'last');

text = text(fi:li);

end
