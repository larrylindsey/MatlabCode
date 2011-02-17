function processJensXLS(filename, outfile, delimiter)

if nargin < 3
    delimiter = sprintf('\t');
else
    delimiter = sprintf(delimiter);
end

res = 1000;

[synstruct conarea dName] = getAreaVsDistance(filename, delimiter);

r = numel(conarea);
minArea = min(conarea);
maxArea = max(conarea);

startArea = floor(res * minArea) / res;
stopArea = ceil(res * maxArea) / res;

testArea = startArea:(1/res):stopArea;

n = numel(testArea);

for i_s = 1:numel(synstruct)
    synstruct(i_s).distmap = calcDistMap(synstruct(i_s), testArea);
end

distMapFull = mkDistMapFull(synstruct, r, n);

writemap(distMapFull, testArea, conarea, dName, outfile);

end

function distMapFull = mkDistMapFull(synstruct, r, n)

%distMapData = zeros(r, n + 1);
distMapFull = zeros(r, n);

for i_s = 1:numel(synstruct)
%    ns = numel(synstruct(i_s).index);
    distMapFull(synstruct(i_s).index, :) = synstruct(i_s).distmap;
end

end

function bigDistMap = calcDistMap(synstruct, testArea)

dx = synstruct.dx;

selectGood = ~isnan(dx);

dxGood = dx(selectGood);
conareaGood = synstruct.conarea(selectGood);
%conareaGood(logical(isnan(conareaGood))) = -inf;

xx = cumsum(dxGood);

n = numel(testArea);

distMap = nan(numel(conareaGood), n);

for i_t = 1:n
    t = testArea(i_t);
    tSelector = logical(conareaGood > t);
    xxT = xx(tSelector);
    
    %distR = [inf diff(xxT)];
    %distL = [diff(xxT) inf];
    
    %dist = min(distR, distL);
    dist = [diff(xxT) nan];
    
    distMap(tSelector, i_t) = dist;
    distMap(not(tSelector), i_t) = nan;   
    
end

bigDistMap = nan(numel(synstruct.conarea), n);
bigDistMap(selectGood, :) = distMap;

end

function [synstruct conarea dendriteName] =...
    getAreaVsDistance(filename, delimiter)

dxTag = 'linNN';
conTag = 'CONarea w/o sym';
nameTag = 'D';

fid = fopen(filename, 'r');

ll = fgetl(fid);
ncell = numel(find(ll == delimiter)) + 1;
controlString = repmat('%s', [1 ncell]);

fclose(fid);

fid = fopen(filename, 'r');

cells = textscan(fid, controlString, 'delimiter', delimiter, ...
    'whitespace', '');

tags = cell(1, numel(cells));

for ii = 1:numel(cells)
    tags{ii} = cells{ii}{1};
end

i_dx = strmatch(dxTag, tags, 'exact');
i_dx = i_dx(1);

i_con = strmatch(conTag, tags, 'exact');
i_con = i_con(1);

i_name = strmatch(nameTag, tags, 'exact');
i_name = i_name(1);

dendriteName = cells{i_name};
dendriteName = {dendriteName{2:end}};

[uniqueDN idn jdn] = unique(dendriteName);

dxCells = cells{i_dx};
conCells = cells{i_con};

dx = zeros(1, numel(dxCells) - 1);
conarea = zeros(1, numel(conCells) - 1);

for ii = 1:(numel(dxCells) - 1)
    dxNum = str2double(dxCells{ii + 1});
    conNum = str2double(conCells{ii + 1});
    
    if numel(dxNum) == 0
        dxNum = nan;
    end
    
    if numel(conNum) == 0
        conNum = nan;
    end
    
    dx(ii) = dxNum;
    conarea(ii) = conNum;
end

fclose(fid);

synstruct = repmat(struct, [1 numel(uniqueDN)]);

for i_s = 1:numel(synstruct)
    selector = find(jdn == i_s);
    synstruct(i_s).name = uniqueDN{i_s};
    synstruct(i_s).dx = dx(selector);
    synstruct(i_s).conarea = conarea(selector);
    synstruct(i_s).index = selector;
end

end

function writemap(distmap, testArea, conArea, dName, outfile)
fid = fopen(outfile, 'w');
n = size(distmap, 2);
%Create the column head cells
topstr = 'CONarea,D,';
comma = ',';
for ii = 1:n
    if ii == n
        comma = '';
    end
    nstr = rpnum2str(testArea(ii), 3);
    topstr = sprintf('%s%s%s', topstr, nstr, comma);
end
fprintf(fid, topstr);
fprintf(fid, '\n');

for i_r = 1:size(distmap, 1)
    comma = ',';
    fprintf(fid, '%s,', rpnum2str(conArea(i_r), 3));
    fprintf(fid, '%s,', dName{i_r});
    for i_n = 1:n
        if i_n == n
            comma = '';
        end
        nstr = rpnum2str(distmap(i_r, i_n), 3);
        fprintf(fid, '%s%s', nstr, comma);
    end
    fprintf(fid, '\n');
    
end

fclose(fid);
end

