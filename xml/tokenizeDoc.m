function [phraseToks, toktype] = tokenizeDoc(instring, delimiters, pairs)

if nargin < 3
    pairs = [''''''; '""'; '()'];
end

if nargin < 2
    delimiters = sprintf(' \t\n');
end

pair1 = pairs(:,1)';
pair2 = pairs(:,2)';
% pair1 = '"''(';
% pair2 = '"'')';

[po pc] = findPairs(instring, pair1, pair2);

pp = zeros(1, numel(po) + numel(pc));
pp(1:2:end) = po;
pp(2:2:end) = pc + 1;
pe = [pp - 1, numel(instring)];
pp = [1 pp];

chunks = {};

for ii = 1:numel(pp)
    sel = pp(ii):pe(ii);
    sel = sel(logical(sel > 0));
    sel = sel(logical(sel <= numel(instring)));
    chunks = cat(1, chunks, {instring(sel)});
end

% Tokenizable chunks have odd indices
% Chunks that should remain intact have even indices.

tokStruct = struct('toks', {});
tokStruct = repmat(tokStruct, size(chunks));
%delimiters = sprintf(' \t\n=');
even = false;

for ii = 1:numel(chunks)

    if even        
        tokStruct(ii).toks.text = chunks{ii};
        tokStruct(ii).toks.type = 1;
    else
        [p tok r] = catchTok(chunks{ii}, delimiters);
        toks = doCat([], p, tok);

        while ~isempty(r)
            [p tok r] = catchTok(r, delimiters);
            toks = doCat(toks, p, tok);
        end        
       
        tokStruct(ii).toks = toks;
    end
    
    even = not(even);
end

tokStruct = cat(1, tokStruct.toks);

for i_t = (numel(tokStruct) - 1):-1:1
    if tokStruct(i_t).type == tokStruct(i_t + 1).type
        tokStruct(i_t).text = [tokStruct(i_t).text tokStruct(i_t + 1).text];
        sel = true(size(tokStruct));
        sel(i_t + 1) = false;
        tokStruct = tokStruct(sel);        
    end
end

phraseToks = {tokStruct.text};
toktype = [tokStruct.type];

end

function toks = doCat(toks, p, tok)
pstr.text = p;
pstr.type = 0;
tstr.text = tok;
tstr.type = 1;

if isempty(pstr.text)
    pstr = [];
end

if isempty(tstr.text)
    tstr = [];
end

toks = cat(1, toks, pstr, tstr);
end

function [po pc] = findPairs(str, pairFH, pairSH)

po = [];
pc = 0;
[pos type] = findFirstStr(str, pairFH);
while ~isempty(pos) && isfinite(pos) && numel(str) > 0
    pos = pos(1);
    po = [po, pos + pc(end)];
    matchChar = pairSH(type);
    str = str((pos + 1):end);
    pos = strfind(str, matchChar);
    if ~isempty(pos)
        pos = pos(1);
        pc = [pc, pos + po(end)];
        str = str((pos + 1):end);
        [pos type] = findFirstStr(str, pairFH);
    end
end

pc = pc(2:end);

end

function [pos type] = findFirstStr(str, list)
pos = zeros(size(list));
for i_t = 1:numel(list)
    p = strfind(str, list(i_t));
    if isempty(p)
        p = inf;
    else
        p = p(1);
    end
    pos(i_t) = p;
end

[pos type] = min(pos);

end


function [p t r] = catchTok(string, delimiters)
[t r] = strtok(string, delimiters);
if isempty(t)
    p = string;
else
    ip = strfind(string, t);
    p = string(1:(ip-1));
end
end

