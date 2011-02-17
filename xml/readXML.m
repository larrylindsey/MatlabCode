function [docp doc] = readXML(infile, control)
% [docp doc] = readXML(infile, [control])
% Reads and parses an XML file.
% infile - a string path to the XML file to read
% control - a control struct containing the following fields:
%   doClean: boolean singleton
%       If set to true, the returned document struct will be stripped of
%       internal control fields and empty elements.
%       Defaults to true.
%   openCloseMatch: boolean array
%       An array of three boolean values.
%       This function considers there to be three types of non-element
%       structures:
%
%       Declaration structures, of the type <?xml ... ?>
%       Meta structures, of the type <! ... >, and
%       Comment structures, of the type <!-- ... -->
%
%       These are generally not parsed beyond simply capturing the text
%       that they enclose.
%       
%       If the corresponding value in this array is set to true, the
%       terminator of this structure will be the '>' matching the first
%       '<'.  The default behavior is for this array to be set to 
%       [false true false], indicating that declaration and comment
%       structure types are terminated by '?>', and '-->', regardless of
%       the text inside of those tags.  Meta tags, however, have no such
%       key characters, and are terminated by a matched '>' character.
%
% docp - a structure representing a parsed XML document.  The structure, if
%   cleaned, will have the following feilds:
%   name: indicates the name of this structure
%   type: indicates the type of this structure
%       'declaration' - declaration structure
%       'meta' - meta structure
%       'comment' - comment
%       'element' - an XML element
%       'text' - a structure that has no name, and contains only textual
%           information in the content field
%       'assignment' - corresponds to an assignment of the type
%           color="red".  The name in this case would be 'color', and the
%           content field would contain '"red"'.
%   pre: contains a cell array of extraneous characters not belonging to
%       any element, but occured before the declaration of the current
%       structure.  For instance, in the case of
%           <picture>
%               <color>red</color>
%           </picture>
%       the struct representing the color element will contain at least a
%       tab and carriage return character in the pre field.
%   attribute: structures defined inside the start-tag of the element in
%       question.
%   content: structures defined in between start- and end-tags
%       of the element in question.
% 
% An example case is included in the comments of this file
%
% doc - an intermediate structure that may be useful for debugging,
%   containing only elements, meta data, comments, and declarations.

%
% Example:
%
% <?xml version="1.0" encoding='UTF-8'?>
% <!-- This is a comment
% -->
% <painting>
%   <img src="madonna.jpg" alt='Foligno Madonna, by Raphael'/>
%   <caption>This is Raphael's "Foligno" Madonna, painted in
%   <date>1511</date>-<date>1512</date>.</caption>
% </painting>
%
% The text above will be converted to a struct array, for instance
% doc = 
% 3x1 struct array with fields:
%     name
%     type
%     pre
%     attribute
%     content
%
% doc(1)
% ans = 
%          name: '?'
%          type: 'declaration'
%           pre: {}
%     attribute: []
%       content: [3x1 struct]
%
% doc(1).content(1)
% ans = 
%          name: ''
%          type: 'text'
%           pre: []
%     attribute: []
%       content: 'xml'
%  doc(1).content(2)
% ans = 
%          name: 'version'
%          type: 'assignment'
%           pre: ' '
%     attribute: []
%       content: '"1.0"'
%  doc(1).content(3)
% ans = 
%          name: 'encoding'
%          type: 'assignment'
%           pre: ' '
%     attribute: []
%       content: ''UTF-8''
%
% doc(2)
% ans = 
%          name: '--'
%          type: 'comment'
%           pre: {}
%     attribute: []
%       content: [1x1 struct]
% doc(2).content
% ans = 
%          name: ''
%          type: 'text'
%           pre: []
%     attribute: []
%       content: 'This is a comment'
%
% doc(3)
% ans = 
%          name: 'painting'
%          type: 'element'
%           pre: {}
%     attribute: []
%       content: [2x1 struct]
% doc(3).content(1)
% ans = 
%          name: 'img'
%          type: 'element'
%           pre: {}
%     attribute: [2x1 struct]
%       content: []
% doc(3).content(2)
% ans = 
%          name: 'caption'
%          type: 'element'
%           pre: {}
%     attribute: []
%       content: [5x1 struct]
%
% and so on.

fid = fopen(infile, 'r');

if fid < 1
    error('Error while opening file %s for reading', infile);
end

if nargin < 2
    control.doClean = true;
    control.openCloseMatch = [false true false];
end
control.fid = fid;

%Read the first structure.
doc = doNextTag(control);

%As long as there are more structures to be read, keep reading them.
while doc(end).continue
   line = doc(end).line;
   doc = cat(1, doc, doNextTag(control, line));
end

fclose(doc(end).control.fid);

%Cleanup.
if control.doClean
    doc = cleanup(doc);
end

%Parse.
docp = parseXML(doc);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstr = doNextTag(control, line)
% control - the control struct
% line - initialize the current line to something.  defaults to empty
if nargin < 2
    line = '';
end

%Init the struct, otherwise there might be cat'ing problems down the line.
fstr.control = control; %Control elements
fstr.tagOpen = [];
fstr.line = line;
fstr.pos = 1;
fstr.continue = true;
fstr.currField = 'attribute';
fstr.state = 0;

fstr.name = ''; %Return elements
fstr.type = '';
fstr.pre = {};
fstr.attribute = {};
fstr.content = {};

%Find the next '<' control character
fstr = nextStartControl(fstr);
%Figure out what type of structure we're dealing with.
fstr = getTagType(fstr);
%Then go and process it.
fstr = processTag(fstr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstr = nextStartControl(fstr)

if isempty(fstr.line)
    fstr = getNextLine(fstr);
end

pos = strfind(fstr.line, '<');

%If we can't find a <, keep going.
while isempty(pos) && fstr.continue
    %if fstr.continue is false, we probably have an error condition.
    
    %Any characters found before the control character, add to the pre
    %field.
    fstr = doCat(fstr, 'pre');
    fstr = getNextLine(fstr);
    pos = strfind(fstr.line, '<');
end

fstr.pos = pos + 1;

fstr = advanceLine(fstr);
%When we've found our control character '<', fstr.pos should indicate the
%first character after it.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = getTagType(fstr)
% I tried to write this to be self-explanatory
% This works by looking at the first non white-space character after a 
% '<' character, which should have been found previously by
% nextStartControl.
%
% look for:
% ? - declaration
% ! -
%     !--, comment
%     otherwise, meta
% otherwise, element.
%
% Technically, having whitespace after a '<' should result in an error
% condition, I think.

line = fstr.line;

nchr = strtok(line);
while isempty(nchr)
    fstr = getNextLine(fstr);
    nchr = strtok(fstr.line);
end

switch nchr(1)
    case '?'
        type = 1;
        fstr.pos = fstr.pos + 1;
        fstr.type = 'declaration';
    case '!'
        if numel(nchr) > 2 && strcmp(nchr(2:3), '--')
            type = 3;
            fstr.pos = fstr.pos + 3;
            fstr.type = 'comment';
        else
            type = 2;
            fstr.pos = fstr.pos + 1;
            fstr.type = 'meta';
        end
    otherwise
        type = 0;
        fstr.type = 'element';
end

tagOpen.type = type;

if type > 0
    tagOpen.char = nchr(1);
else
    tagOpen.char = '';
end

fstr.tagOpen = tagOpen;

fstr = advanceLine(fstr);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = processTag(fstr)
%Dispatch to the correct function to handle our structure type.
switch fstr.tagOpen.type
    case 0
        fstr = parseTagInternal(fstr);
    case 1
        fstr = parseComment(fstr, '?', fstr.tagOpen.type);
    case 2
        fstr = parseComment(fstr, '', fstr.tagOpen.type);
    case 3
        fstr = parseComment(fstr, '--', fstr.tagOpen.type);
    otherwise
        error('Unknown tag type %d', fstr.tagOpen.type);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstr = parseComment(fstr, term, type)
%Grabs all of the text in a non-data-y tag.
% term - some characters to look for before the '>', to determine when the
%    structure is termindated.
% type - an integer indicating the structure type.
fstr.name = term;

%If we need to match the close '>' to the open '<':
if fstr.control.openCloseMatch(type)
    % This is used especially in the case of meta declarations, which
    % includes DTD's.
    % You might see something like
    % <!bla [ <! more bla >
    %   <!foo >
    %   <!then bar> ] >
    % 
    % These tags are not terminated in the regular way, so we just look for
    % the '>' that matches the first '<', by counting:
    % <  < < >  < > >  >
    % 1  2 3 2  3 2 1  0<-end    
    
    count = 1;
    fstr.currField = 'content';
    
    fstr = advanceLine(fstr);
    
    while count > 0
        posClose = strfind(fstr.line, '>');
        posOpen = strfind(fstr.line, '<');

        if isempty(posClose)
            posClose = inf;
        end
        if isempty(posOpen)
            posOpen = inf;
        end

        [pos type] = min([posClose posOpen]);

        if isinf(pos)
            fstr = doCat(fstr, 'content');
            fstr = getNextLine(fstr);
        else
            fstr = doCat(fstr, 'content', fstr.line(1:pos));
            fstr.pos = pos + 1;
            fstr = advanceLine(fstr);
            if type == 1
                count = count - 1;
            else
                count = count + 1;
            end
        end
    end
    
    closeGrel = findstr(fstr.content{end}, '>');
    closeGrel = closeGrel(end);
    fstr.content{end} = fstr.content{end}(1:(closeGrel - 1));
else
    % When not matching by counting opens/closes, we specifically look for
    % a certain character string.  In the case of a declaraction, this
    % should be '?>', and in the case of a comment, it is '-->'.  Any
    % characters that occur before these strings are stored in
    % fstr.content, but otherwise ignored.
    
    term = [term '>'];
    fstr.currField = 'content';
    fstr = advanceLine(fstr);
    fstr = doCat(fstr);
    pos = strfind(fstr.line, term);
    while isempty(pos)
        fstr = getNextLine(fstr);
        fstr = doCat(fstr, 'content');
        pos = strfind(fstr.line, term);
    end

    lastline = fstr.content{end};
    lastline = lastline(1:(pos - 1));
    fstr.content{end} = lastline;
    fstr.pos = pos + numel(term);
    fstr = advanceLine(fstr);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstr = parseTagInternal(fstr)
tok = strtok(fstr.line);
while isempty(tok)
    fstr = getNextLine(fstr);
    tok = strtok(fstr.line);
end

% Chicken and egg!
% We must first find the name of this tag before we can really do any
% processing, but we have to reject any transitional characters that may
% appear in the same token.
grpos = strfind(tok, '>');
slpos = strfind(tok, '/');

if numel(grpos > 1)
    grpos = grpos(1);
end
if numel(slpos > 1)
    slpos = slpos(1);
end

% Find the first > or / in the token
% There will not be one in the case of a tag like
% <tagname foo=bar/>
% but we might expect one for
% <tagname/>, or
% <tagname> foo=bar </tagname>
%
% We want just the "tagname" part.
tokl = min([grpos,slpos]);
if isempty(tokl)
    tagName = tok(1:end);
else
    tagName = tok(1:(tokl - 1));    
end

fstr.name = tagName;

fstr.pos = strfind(fstr.line, tagName) + numel(tagName);

% OK. Now we have the tag name, and its position.  Time to remove it from
% the line we're looking at.
fstr = advanceLine(fstr);

% Now we begin processing in earnest.
% We do the actual parsing later.
% Processing here takes the form of finding a transition and handling it,
% until there are no more transitions to handle.
trans = findTransition(fstr.line);

fstr = handleTrans(fstr, trans);

while ~trans.terminator
    trans = findTransition(fstr.line);
    fstr = handleTrans(fstr, trans);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = handleTrans(fstr, trans)
%If trans.type is zero, we just concatenate the current line, then move on.
% In each other case, we concatenate up until the transition, then continue
% to process it.
switch trans.type    
    case 0 
        % No terminator.
        fstr = doCat(fstr);
        fstr = getNextLine(fstr);
    case 1
        % Found >
        % This should mean that we're transitioning out of an xml start-tag
        % with a corresponding end-tag.
        if fstr.state == 0
            line = fstr.line;
            line = line(1:(trans.pos - 1));
            fstr = doCat(fstr, [], line);
            fstr.pos = trans.pos + 1;
            fstr = advanceLine(fstr);
            fstr.state = 1;
            fstr.currField = 'content';
        else
            error('Unexpected tag state');
        end
    case 2
        % Found <, indicating the beginning of a start-tag
        line = fstr.line;
        line = line(1:(trans.pos - 1));
        fstr = doCat(fstr, [], line);
        substr = doNextTag(fstr.control, fstr.line);
        fstr = doCat(fstr, [], substr);
        fstr.line = substr.line;
        fstr.pos = substr.pos;
    case 3        
        % Found />, the end of a start-tag that terminates the element.
        line = fstr.line;
        line = line(1:(trans.pos - 1));
        fstr = doCat(fstr, [], line);
        fstr.pos = trans.pos + 2;
        fstr = advanceLine(fstr);
    case 4
        % Found </, the beginning of an end-tag.
        line = fstr.line;
        line = line(1:(trans.pos - 1));
        fstr = doCat(fstr, [], line);
        %We expect the end-tag to match the start-tag.  We want to find the
        %end of the end-tag before we can move on.
        
        %First try to find the name of the tag.
        pos = strfind(fstr.line, fstr.name);        
        cnt = 0;
        while isempty(pos)
            fstr = getNextLine(fstr);
            pos = strfind(fstr.line, fstr.name);
            cnt = cnt + 1;
            %Try twice, then give up.
            if cnt > 1
                error('Could not find terminator tag');
            end            
        end
        pos = pos(1);
        fstr.pos = pos + numel(fstr.name);
        fstr = advanceLine(fstr);
        %Now, we're past the tag name, try to find the terminating >
        pos = strfind(fstr.line, '>');
        cnt = 0;
        while isempty(pos)
            fstr = getNextLine(fstr);
            pos = strfind(fstr.line, '>');
            cnt = cnt + 1;
            %Again, try twice, then give up.
            if cnt > 1
                error('Missing control character: >');
            end
        end
        pos = pos(1);
        fstr.pos = pos + 1;
        fstr = advanceLine(fstr);
        %Now, the structure has terminated, and we can move on. The
        %position should now be just one index beyond the end of the
        %end-tag we just found.
    otherwise
        error('Unknown transition type %d', trans.type);
end

if isempty(strtok(fstr.line))
    fstr = getNextLine(fstr);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trans = findTransition(line)
% Tag transitions are defined by > or < characters.  Slashes indicate the
% termination of a structure.
gr = strfind(line, '>');
lt = strfind(line, '<');
grterm = strfind(line, '/>');
ltterm = strfind(line, '</');

if isempty(gr); gr = inf; end
if isempty(lt); lt = inf; end
if isempty(grterm); grterm = inf; end
if isempty(ltterm); ltterm = inf; end

gr = gr(1);
lt = lt(1);
grterm = grterm(1);
ltterm = ltterm(1);

if ltterm <= lt
    lt = inf;
end

% 1 - ends tag open
% 2 - beginning of a subtag
% 3 - ends tag open, terminates current tag
% 4 - beginning of a tag terminator
[trans.pos trans.type] = min([gr lt grterm ltterm]);

% 0 - no transitions found on this line
if isinf(trans.pos)
    trans.type = 0;
end

trans.terminator = (trans.type >=3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = getNextLine(fstr)
%fstr should contain a control struct that has, among others, a field
%containing the fid for the file we're reading.
line = fgetl(fstr.control.fid);
while isempty(strtok(line))
    line = fgetl(fstr.control.fid);
end

fstr.line = line;
fstr.pos = 1;

% If the file is done, fgetl will return a noncharacter, -1 to be specific.
% In this case, we have to indicate that processing has finished.  This is
% done by setting fstr.continue to false.
if ~ischar(line)
    fstr.continue = false;
%    fclose(fstr.control.fid);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = advanceLine(fstr)
% This function takes care of the bookkeeping we need to do to advance the
% position of the line.  fstr.pos should be set to the position of the
% character that should begin the new line.  In other words, if we want
% to set fstr.line to fstr.line(3:end), then fstr.pos should be 3.
%
% If the resulting string has no tokens in it, we make a call to
% getNextLine until we find one that does.
%
% TODO: figure out what to do with non whitespace characters.
if numel(fstr.pos) > 1
    fstr.pos = fstr.pos(1);
end
fstr.line = fstr.line(fstr.pos:end);
if isempty(strtok(fstr.line))
    fstr = getNextLine(fstr);
end
fstr.pos = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = doCat(fstr, field, text)
% concatenate the contents of a given field, by adding a new element in a
% cell array.
%

if nargin < 3
    text = fstr.line;
end

if nargin < 2 || isempty(field)
    field = fstr.currField;
end

if ~isempty(text)
    fstr.(field) = cat(1, fstr.(field), {text});
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fstr = cleanup(fstr)
% Removes control fields and empty content from the struct tree.
%
delFields = {'control', ...
    'continue', ...
    'line', ...
    'tagOpen', ...
    'state', ...
    'currField', ...
    'pos'};

checkFields = {'content', 'attribute'};

for i_s = 1:numel(fstr)
    for i_f = 1:numel(checkFields)
        field = fstr(i_s).(checkFields{i_f});
        keep = true(size(field));
        for i_c = 1:numel(field)
            if isstruct(field{i_c})
                field{i_c} = cleanup(field{i_c});
            elseif ischar(field{i_c})
                if isempty(strtok(field{i_c}))
                    keep(i_c) = false;
                end
            end
        end
        field = {field{keep}};
        fstr(i_s).(checkFields{i_f}) = field;
    end
end

fstr = rmfield(fstr, delFields);
end