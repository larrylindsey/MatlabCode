function doc = parseXML(doc)
% doc = parseXML(doc)
%
% Parses a struct tree of the type returned in the second output of readXML
%

% doc should be a struct array.  Run parseDoc on each element.
for i_s = 1:numel(doc)
    doc(i_s) = parseDoc(doc(i_s));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = parseDoc(fstr)
% Create the template.
template = fstr(1);
fn = fieldnames(template);
for i_f = 1:numel(fn);
    field = fn{i_f};
    template.(field) = [];
end

%Data fields.
fields = {'content', 'attribute'};

% If we get a text or assignment structure type, skip it.
% NOTE: Upon re-reading this code, I'm not sure that this conditional is
% necessary, but I'm leaving it in just in case.
if ~strcmp(fstr.type, 'text') && ~strcmp(fstr.type, 'assignment')
    for i_f = 1:numel(fields)
        %Parse the contents of the content field, then the attribute field.
        %There is no significance to the order I did this in.
        fstr.(fields{i_f}) = doParse(fstr.(fields{i_f}), template);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newCArray = doParse(cArray, template)

% Data should be passed in as a cell array.
% The template variable gives us a struct to push the parsed data into.

newCArray = {};

i_c = 1;

while i_c <= numel(cArray)
    %Get the current element
    el = cArray{i_c};
    if isstruct(el)
        %If it is a struct, then it should have the same structure as the
        %doc struct array passed in. Parse it as such.
        el = parseDoc(el);
        %cat the result into the output cell array
        newCArray = cat(1, newCArray, {el});
        i_c = i_c + 1;
    elseif ischar(el)
        %If the element is a string, concatenate it with all of the
        %subsequent strings in the cell array, then parse it as a string.
        parseStr = '';
        while ischar(el) && i_c <= numel(cArray)
            parseStr = sprintf('%s\n%s', parseStr, el);
            i_c = i_c + 1;
            if i_c <= numel(cArray)
                el = cArray{i_c};
            end
        end
        strParsed = doParseString(parseStr, template);
        newCArray = cat(1, newCArray, strParsed);
    else
        error('Unexpected data type in attribute or conent field');
    end
    
end

% Our new cell array should contain structs with identical fieldnames, so
% we may cat them.  This wasn't done from the beginning for debuggability
% purposes.
newCArray = cat(1,newCArray{:});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parseArray = doParseString(instring, template)
%Parse by an '=' as well as by whitespace.

[toks types] = tokenizeDoc(instring, sprintf(' \t\n='));

isToken = logical(types);

% toks comes back as an alternating array of tokens and token separators
% (one of ' \t\n=').  We'd like the first value in the tok array to be an
% actual text-carrying token, then toks{1:2:end} should correspond to text
% tokens, and toks{2:2:end} should correspond to separators.
while ~isempty(isToken) && not(isToken(1))
    toks = {toks{2:end}};
    isToken = isToken(2:end);
end

%Collect all of the assignments

if isempty(isToken)
    % If there are no tokens, return an empty struct.
    parseArray = reshape(template, [1, 0]);
else
    % Initialize the first parsed array struct to a text XML structure.
    parseArray = template;
    parseArray.type = 'text';
    parseArray.name = '';
    for i_t = 1:2:numel(toks)
        
        parseArray(end).content = toks{i_t};
        if i_t > 1
            % if i_t is 1, then there is no inter-token to record.
            if strcmp(parseArray(end).type, 'text')                
                parseArray(end).pre = toks{i_t - 1};
            else
                %Assume type is 'assignment'
                parseArray(end).attribute = toks{i_t - 1};
            end
        end
        if numel(toks) - i_t >= 2
            if strfind(toks{i_t + 1}, '=')
                % if the current inter-token contains an equal sign, this
                % should be an assignment.
                % The name of this struct gets the current token.
                % The content contains the token too, right now, but it
                % should be set correctly in the next iteration.
                parseArray(end).name = toks{i_t};
                parseArray(end).type = 'assignment';
            else
                % If the inter-token did not contain an equal sign, leave
                % the current struct as-is, then cat a new one to the end
                % of the array.
                parseArray = cat(1, parseArray, template);
                parseArray(end).type = 'text';
                parseArray(end).name = '';
            end
        end
    end
end

%Now, cat all of the consecutive text-only strings together.
keepsel = true(size(parseArray));
for i_t = (numel(parseArray) - 1):-1:1
    if strcmp(parseArray(i_t).type, 'text') && ...
            strcmp(parseArray(i_t + 1).type, 'text')
        parseArray(i_t).content = [parseArray(i_t).content ...
            parseArray(i_t + 1).pre parseArray(i_t + 1).content];
        keepsel(i_t + 1) = false;
    end        
end

parseArray = parseArray(keepsel);
end