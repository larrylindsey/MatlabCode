function xmlstr = generateXML(in, name, lvl)
%xmlstr = generateXML(in)
%
%  Output the variable "in" as an XML string; intended for structures. When
%  cell arrays or sub-structures are encountered, a foreign key is added in
%  the parent object. The child cell array or structure recive a
%  corresponding primary key. This allows importing the xml in database
%  systems that generates new tables for each cell array and substructure
%  while keeping the relationships...
%  The datatypes supported are: char, numeric, struct, cell, logical/boolean
%
%  This XML format generates a HUGE overhead. If you do not intend to import
%  the data in a database, but only to save in human-readable form the data
%  for further use in MatLab, it is advised to use the more efficient XML
%  Toolbox by Marc Molinari, available on Matlab exchange website.
%  generateXML is a one-way function to export structures, it does not
%  feature an import fonctionality [yet]... to reimport data, use SQL and
%  the CSE SQL Library, available on Matlab exchange website.
%
% Written by L.Cavin, 2004

global uniqueID;

% check input parameters
if (nargin<1)
  error('Sorry, but I need something to export...');
end
if nargin<2
    name = inputname(1);
    if isempty(name),  name = 'root'; end
end
if (nargin<3) || isempty(lvl), lvl = 0; end

if exist('uniqueID') && ~isempty(uniqueID) && lvl == 0
    % do the things cleanly, as we use the global workspace:
    warning('CSE:generateXML', 'Global Variable "uniqueID" will be erased by generateXML.');
end
if ~exist('uniqueID') || isempty(uniqueID)
    uniqueID = 0;
end

xmlstr = '';
depth = repmat(' ', 1, 3*(lvl+1));
new_line = sprintf('\n');

if lvl == 0
    xmlstr = ['<xml>' new_line];
    xmlstr = [xmlstr '<!-- source="MatLab" encoded="generateXML by CSE, 2004" -->' new_line];
    xmlstr = [xmlstr '<!-- generated ' datestr(now) ' -->' new_line];
end

switch lower(class(in))
    case 'double'
        if isreal(in)
            for cntr = 1:prod(size(in))
                % write wrapper
                nme = [name, '_' num2str(cntr)];
                xmlstr = [xmlstr, depth, '<', nme, ' size="' deblank(sprintf('%d ', size(in))) '" type="' class(in) '">'];
                content = sprintf('%0.16g ', in(cntr));
                xmlstr = [xmlstr, content(1:end-1), '</', nme, '>', new_line];
            end
        else
            % write wrapper
            str = feval('generateXML', real(in), [name '_real'], lvl);
            xmlstr = [xmlstr, sprintf('%s', str)];
            str = feval('generateXML', imag(in), [name '_imaginary'], lvl);
            xmlstr = [xmlstr, sprintf('%s', str)];
        end
    case 'boolean'
        for cntr = 1:prod(size(in))
            % write wrapper
            nme = [name, '_' num2str(cntr)];
            xmlstr = [xmlstr, depth, '<', nme, ' size="' deblank(sprintf('%d ', size(in))) '" type="' class(in) '">'];
            content = sprintf('%0.16g ', in(cntr));
            xmlstr = [xmlstr, content(1:end-1), '</', nme, '>', new_line];
        end
    case 'char'
        % write wrapper
        xmlstr = [xmlstr, depth, '<', name, ' size="' deblank(sprintf('%d ', size(in))) '" type="' class(in) '">'];
        % substitute XML reserved characters with their html equivalent
        content = htmlize(in(:)');
        xmlstr = [xmlstr, content, '</', name, '>', new_line];
    case 'struct'	
        % write ID-link with new table
        if lvl > 0
            blp = uniqueID;
            uniqueID = uniqueID + 1;
            xmlstr = [xmlstr, depth, '<link_' name ' type="Foreign Key" table="' name '">' int2str(blp) '</link_' name '>' new_line];
        end
        flds = fieldnames(in);
        for cntr = 1:prod(size(in))
            % write wrapper
            xmlstr = [xmlstr, depth, '<', name, ' index="' num2str(cntr) '" size="' deblank(sprintf('%d ', size(in))) '" type="' class(in) '">' new_line];
            if lvl > 0
                xmlstr = [xmlstr, depth, '   <' name 'ID type="Primary Key">' int2str(blp) '</' name 'ID>' new_line];
            end
            for n = 1:length(flds)
                % get content
                content = getfield( in(cntr), flds{n} );
                % write content
                str = feval('generateXML', content, flds{n}, lvl+1);
                xmlstr = [xmlstr, sprintf('%s', str)];
            end
            % finishwrapper
            xmlstr = [xmlstr, depth(1:end) '</', name, '>' new_line];
        end 
    case 'cell'
        % write ID-link with new table
        blp = uniqueID;
        uniqueID = uniqueID + 1;
        xmlstr = [xmlstr, depth, '<link_' name ' type="Foreign Key" table="' name '">' int2str(blp) '</link_' name '>' new_line];
        % write wrapper
        xmlstr = [xmlstr, depth, '<', name, ' size="' deblank(sprintf('%d ', size(in))) '" type="' class(in) '">' new_line];
        xmlstr = [xmlstr, depth, '   <' name 'ID type="Primary Key">' int2str(blp) '</' name 'ID>' new_line];
        for n = 1:length(in);
            content = in{n};
            % write content
            xmlstr = [xmlstr, feval('generateXML', content, [name '_' int2str(n)], lvl+1)];
        end
        % finishwrapper
        xmlstr = [xmlstr, depth(1:end-2), '  </', name, '>' new_line];
    otherwise
        disp(['Current Type: ', att.type]);
        error(['Use only implemented types: double (also complex), char, struct, cell, boolean.']);
end

if lvl == 0
    xmlstr = [xmlstr '</xml>'];
    clear global uniqueID;
end
    
function in = htmlize(in)
% regular expressions would be more suited, in the form:
%   in = regexprep(in, '&', '&amp;');
%   in = regexprep(in, '<', '&lt;');
%   in = regexprep(in, '>', '&gt;');
%   in = regexprep(in, '''', '&apos;');
%   in = regexprep(in, '"', '&quot;');
% But as a french-speaking living in a german speaking region, I must
% support the use of יאטח... as well as üצה... so I cannot use this
% function and hence this code blows up from 5 lines to 25 lines...
% (REGEXPREP does not support international character sets.)
stridx = findstr(in, '&');
for cntr=1:length(stridx)
  in = [in(1:stridx(cntr)-1), '&amp;', in(stridx(cntr)+1:end)];
  stridx = stridx+4;
end
stridx = findstr(in, '<');
for cntr=1:length(stridx)
  in = [in(1:stridx(cntr)-1), '&lt;', in(stridx(cntr)+1:end)];
  stridx = stridx+3;
end
stridx = findstr(in, '>');
for cntr=1:length(stridx)
  in = [in(1:stridx(cntr)-1), '&gt;', in(stridx(cntr)+1:end)];
  stridx = stridx+3;
end
stridx = findstr(in, '''');
for cntr=1:length(stridx)
  in = [in(1:stridx(cntr)-1), '&apos;', in(stridx(cntr)+1:end)];
  stridx = stridx+5;
end
stridx = findstr(in, '"');
for cntr=1:length(stridx)
  in = [in(1:stridx(cntr)-1), '&quot;', in(stridx(cntr)+1:end)];
  stridx = stridx+5;
end
