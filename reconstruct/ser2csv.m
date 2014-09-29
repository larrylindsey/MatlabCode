function ser2csv(secdoc, rexp, fname)
% ser2csv(secdoc, rexp, fname)
%
%  Do some SER calculations in a Reconstruct project and write the results
%     to a csv file.
%
% secdoc - the sections returned from readReconstruct
% rexp - a regular expression identifying the SER traces in question
% fname - the file to write out

 fprintf('enumerating contours\n');
 
 names = enumerateContours(secdoc);
 sernames = strsearch(names, rexp);

  n = max([secdoc.index]);

 
 ser = repmat(struct, size(sernames));
 
 fprintf('calculating areas, counting cross sections\n');
 
 for ii = 1:numel(sernames)
     ser(ii).name = sernames{ii};
     ser(ii).contour = extractContour(secdoc, ser(ii).name);
     [c a] = countContourXSec(ser(ii).contour);
     a = cat(2, a, zeros(size(a, 1), n - size(a, 2)));
     c = cat(2, c, zeros(1, n - size(c, 2)));
     ser(ii).area = a;
     ser(ii).ct = c;
 end


fprintf('writing\n');

if nargin > 2
    fid = fopen(fname, 'w');
else
    fid = -1;
end

fprintSER(ser, n, fid);

if fid > 0
    fclose(fid);
end

end

function fprintSER(ser, n, fid)

headcstr = '"",';
numblock = zeros(1, n);
numblock(:) = 1:n;
for ii = 1:numel(ser)
    headcstr = cat(2, headcstr, '"', ser(ii).name, '"');
    headcstr = cat(2, headcstr, repmat(',""',[1 1+size(ser(ii).area, 1)]));
    headcstr = cat(2, headcstr, ',');
    numblock = cat(1, numblock, ser(ii).ct, ser(ii).area, ...
        sum(ser(ii).area, 1));
end
dopr(fid, [headcstr '\n']);

numcstr = repmat('%g,', [1 size(numblock, 1)]);
numcstr = cat(2, numcstr, '\n');

for ii = 1:size(numblock, 2)
    numline = numblock(:, ii);
    arg = num2cell(numline);
    dopr(fid, numcstr, arg{:});
end

end

function dopr(fid, cstr, varargin)

if fid > 0
    fprintf(fid, cstr, varargin{:});
else
    fprintf(cstr, varargin{:});
end
end
