function output = writeXML(doc, file, control)
%writeXML(doc, file, [control])
%
% A quick and dirty file-writer for structures of the type returned by
% readXML or parseXML.
%
% doc - the doc struct to write.
% file - the output file name
% control - a struct that controls how the writing is done.
%           call this function with no arguments to return the default
%           struct.
%           

fid = fopen(file, 'w');

if fid < 1
    error('Could not open file %s for writing', file);
end

if nargin < 3

    methods.element = @convertElement;
    methods.meta = @convertMeta;
    methods.declaration = @convertDeclaration;
    methods.comment = @convertComment;
    methods.text = @convertText;
    methods.assignment = @convertAssignment;
    control.methods = methods;
    control.catText = @catText;
    
    if nargin < 1
        output = control;
        return;
    end
end

doctext = dispatch(doc, control);

if nargout > 0
    output = fprintf(fid, '%s', doctext);
else
    fprintf(fid, '%s', doctext);
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = dispatch(str, control, text)

if nargin < 3
    text = '';
end

if ischar(str)
    text = control.catText(text, str);
else
    methodstruct = control.methods;
    methods = fieldnames(methodstruct);
    for i_s = 1:numel(str)
        type = str(i_s).type;
        if isempty(strmatch(type, methods, 'exact'))
            error('No handler for type %s', type);
        end
        text = control.catText(text,...
            methodstruct.(type)(str(i_s), control));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = catText(first, second)

% Condition second

%Make sure a '<' is immediately followed by a something else.
grefind = strfind(second, '<');
for i_g = 1:numel(grefind)
    grefind = grefind(i_g);
    pre = second(1:grefind);
    post = second((grefind + 1):end);
    tok = strtok(post);
    loc = strfind(post, tok);
    post = post(loc:end);
    second = [pre post];
    grefind = strfind(second, '<');    
end

if isempty(first)
    text = second;
else
    %Insert a carriage return every so often.
    if numel(first) + numel(second) > 80
        inter = sprintf('\n');
    else
        inter = ' ';
    end
    text = [first, inter, second];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = convertElement(str, control)
text = [sprintf('\n') '<' str.name];
text = dispatch(str.attribute, control, text);
if isempty(str.content)
    text = [text '/>' sprintf('\n')];
else
    text = [text '>'];
    text = dispatch(str.content, control, text);
    text = control.catText(text, ['</' str.name '>' sprintf('\n')]);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = convertAssignment(str, control)
text = control.catText('',[str.name, str.attribute , str.content]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = convertText(str, control)
text = control.catText('', str.content);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = convertMeta(str, control)
if (numel(str.content) >= 1) && strcmp(str.content(1).type, 'text')
    text = ['<!', str.content(1).content];
else
    text = '<!';
end
text =  dispatch(str.content(2:end), control, text);
text = [text , '>'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = convertDeclaration(str, control)
if (numel(str.content) >= 1) && strcmp(str.content(1).type, 'text')
    text = ['<?', str.content(1).content];
else
    text = '<?';
end
text =  dispatch(str.content(2:end), control, text);
text = [text , '?>'];
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text = convertComment(str, control)
text = control.catText([sprintf('\n') '<!--'],...
    dispatch(str.content, control));
text = [text , '-->'];
end