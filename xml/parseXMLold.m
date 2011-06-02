function [xmlStruct outPos tagName] = parseXML(inXML, inpos, xmlStruct)

ws = ' <>="/';


if nargin < 2
    r = 1;
    c = 1;
else
    r = inpos(1);
    c = inpos(2);
end

if nargin < 3
    xmlStruct = struct;
end

currLine = inXML{r};

if currLine(c) == '<'
    [tagName pos] = nextTok(inXML, r, c+1, ws);

    r = pos(1);
    c = pos(2);
else 
    tagName = '';
end



[xmlStruct outPos] = populateStruct(inXML, r, c, xmlStruct);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xmlStruct outPos] = populateStruct(inXML, r, c, xmlStruct)

if nargin < 4
    xmlStruct = struct;
end

[endPos termHasSlash] = findTerminator(inXML, r, c);
[nextPos nextHasSlash] = findOpen(inXML, r, c);


if comesBefore(endPos, nextPos)
    xmlStruct = placeField(xmlStruct, getTags(inXML, [r c], endPos), ...
        'internal');
    %xmlStruct = placeInternalField(xmlStruct,...
    %    getTags(inXML, [r c], endPos));
    
    if termHasSlash
        outPos = endPos;
    else
        [xmlStruct outPos] = getTags(inXML, endPos, nextPos, xmlStruct);
    end

else
%    xmlStruct = placeInternalField(xmlStruct, ...
%        getTags(inXML, [r c], nextPos));
    xmlStruct = placeField(xmlStruct, getTags(inXML, [r c], nextPos),...
        'internal');
    [subStruct subPos tagName] = parseXML(inXML, nextPos);
    xmlStruct = placeField(xmlStruct, subStruct, tagName);
    [xmlStruct outPos] = parseXML(inXML, subPos, xmlStruct);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos isSlash] = findTerminator(inXML, r, c)
pos = findChar(inXML, r, c, '>');
if c < 2
    isSlash = false;
else
    currLine = inXML{pos(1)};
    isSlash = currLine(pos(2) - 1) == '/';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos isSlash] = findOpen(inXML, r, c)
pos = findChar(inXML, r, c, '<');
currLine = inXML{r};
if pos(2) == numel(currLine)
    isSlash = false;
else
    isSlash = currLine(pos(2) + 1) == '/';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pos = findChar(inXML, r, c, e)
found = false;
while ~found
    currLine = inXML{r};
    currLine = currLine(c:end);
    c = find(currLine == e, 1, 'first');
    if isempty(c)
        c = 1;
        r = r + 1;
    else
        found = true;
    end    
end

pos = [r c];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xmlStruct pos] = getTags(inXML, start, finish, xmlStruct)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xmlStruct = placeFields(xmlStruct, subStruct)

fields = fieldnames(subStruct);

for i_f = 1:numel(fields)
    currField = fields{i};    
    match = strmatch(currField, fieldnames(xmlStruct), 'exact');
    
    if isempty(match)
        datum = {subStruct.(currField)};
    else
        datum = {xmlStruct.(fields{i_f}){:}, subStruct.(currField)};
    end

    xmlStruct.(currField) = datum;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xmlStruct = placeInternalField(xmlStruct, subStruct)

addStruct.internal = subStruct;

xmlStruct = placeFields(xmlStruct, addStruct);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tok pos] = nextTok(xml, r, c, ws)
tok = '';
count = 0;
origC = c;
while isempty(tok)
    currLine = xml{r};
    currLine = currLine(c:end);
    tok = strtok(currLine, ws);
    r = r + 1;
    c = 1;
    count = count + 1;
end

c = numel(tok);

if count == 1
    c = c + origC;
end

pos = [r - 1, c];

end