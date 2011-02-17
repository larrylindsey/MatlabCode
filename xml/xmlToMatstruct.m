function [mstr name] = xmlToMatstruct(xmldoc)

[mstr name] = dispatch(xmldoc(end));
mstr.pre = xmldoc(1:2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mstr name] = dispatch(xmldoc)

mstr = struct;
xml = xmldoc(1);

if strcmp(xml.type, 'element')
    for i_a = 1:numel(xml.attribute)
        [attr name] = dispatch(xml.attribute(i_a));
        mstr = addfield(mstr, name, attr);
    end
    for i_c = 1:numel(xml.content)
        [cont name] = dispatch(xml.content(i_c));
        mstr = addfield(mstr, name, cont);
    end
    name = xml.name;
elseif strcmp(xml.type, 'assignment')
    %mstr = addfield(mstr, xml.name, xml.content);
    name = xml.name;
    mstr = xml.content;
elseif strcmp(xml.type, 'comment')
    %mstr = addfield(mstr, 'comment', xml.content);
    mstr = xml.content;
    name = 'comment';
elseif strcmp(xml.type, 'text')
    %mstr = addfield(mstr, 'text', xml.content);
    mstr = xml.content;
    name = 'text';
else
    error('Unrecognized type %s', xml.type);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mstr =  addfield(mstr, name, content)

%prepare content
if ischar(content)
    content = stripquotes(content);
    numcontent = str2num(content);
    if ~isempty(numcontent)
        content = numcontent;
    end
end

if isfield(mstr, name)
    field = concatenate(mstr.(name), content);
else
    field = content;    
end

mstr.(name) = field;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = stripquotes(str)
str = str(logical(str ~= '"'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = concatenate(field, val)
if ischar(field)
    field = {field};
end

if iscell(field)
    field = {field{:} val};   
elseif isstruct(field)
    field(end + 1) = field(1);
    
    oldfields = fieldnames(field);
    newfields = fieldnames(val);
    
    for i_f = 1:numel(newfields)
        name = newfields{i_f};
        field(end).(name) = val.(name);
    end
    
    nilfields = setdiff(oldfields, newfields);
    
    for i_f = 1:numel(nilfields)
        name = nilfields{i_f};
        field(end).(name) = [];
    end
else
    field(end + 1) = val;
end
end
