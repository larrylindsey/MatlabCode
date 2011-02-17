function outCoords = splitImage(imageFile, rc, ovlp, doWrite, list, z)

if nargin < 5
    list = '';
    if nargin < 4
        doWrite = true;
    end
end

if ~doWrite && logical(nargout == 0)
    return;
end

fid = fopen(list, 'a');


if ischar(imageFile)
    disp('Reading File');
    im = imread(imageFile);
    disp('Done reading');
else
    im = imageFile;
    imageFile = 'output.png';
end

sz = size(im);

rtSz = round(sz(1:2) ./ rc);

olSz = round(rtSz * ovlp);

iDot = find(imageFile == '.', 1, 'last');

if isempty(iDot)
    rootName = imageFile;    
    ext = '.png';
else
    rootName = imageFile(1:(iDot - 1));
    ext = imageFile(iDot:end);
end

rowStartPos = 1:rtSz(1):sz(1);
colStartPos = 1:rtSz(2):sz(2);

if nargout > 0
    outCoords = struct;
    outCoords = repmat(outCoords, [rc(1) rc(2)]);
end

for i_r = 1:rc(1)
    for i_c = 1:rc(2)
%        fprintf('Row %g Column %g\n', i_r, i_c);
        rowStart = max(1, rowStartPos(i_r) - olSz(1));
        rowEnd = min(sz(1), rowStartPos(i_r) + rtSz(1) + olSz(1));
        colStart = max(1, colStartPos(i_c) - olSz(2));
        colEnd = min(sz(2), colStartPos(i_c) + rtSz(2) + olSz(2));
        
        imFrag = im(rowStart:rowEnd, colStart:colEnd, :);
        
        imName = [rootName sprintf('_row%g_col%g', i_r, i_c) ext];

        dirName = sprintf('r%g_c%g', i_r, i_c);
        outName = [dirName '/' imName];

        if doWrite 
            if ~isdir(['./', dirName]);
                mkdir('./', dirName);
            end
    
            imwrite(imFrag, outName);
        end
        
        if fid > 0
            fprintf(fid, '%s\t%g\t%g\t%g\n', outName,...
                colStart, rowStart, z);
        end
        
        if nargout > 0
            outCoords(i_r, i_c).file = outName;
            outCoords(i_r, i_c).r = rowStart;
            outCoords(i_r, i_c).c = colStart;
            outCoords(i_r, i_c).im = imFrag;
        end
        
        clear imFrag;
    end
end

if fid > 0
    fclose(fid);
end

end
