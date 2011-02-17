function displayXML(doc, depth, ntab)

if nargin < 3 || isempty(ntab)
    ntab = 2;
end

if nargin < 2 || isempty(depth)
    depth = inf;
    fprintf('Displaying XML tree using unlimited depth\n');
end

doDisplay(doc, ntab, depth + 1, 1);

end

function doDisplay(doc, ntab, dmax, d)

if d > dmax
    return;
else
    tab = repmat(' ', [1, ntab]);
    sprefix = repmat(tab, [1, d - 1]);
    lprefix = repmat(tab, [1, d]);

    if isstruct(doc)        
        for i_doc = 1:numel(doc)
            if strcmp(doc(i_doc).type, 'text')
                disptext = fixtext(doc(i_doc).content);
                fprintf('%s%s\n', sprefix, disptext);
            elseif strcmp(doc(i_doc).type, 'assignment')
                disptext = fixtext(doc(i_doc).content);
                fprintf('%s%s = %s\n', sprefix, doc(i_doc).name, disptext);
            else
                fprintf('%s%s %s\n', sprefix, doc(i_doc).type, ...
                    doc(i_doc).name);
                if ~isempty(doc(i_doc).attribute) && d + 1 <= dmax
                    fprintf('%sAttributes:\n', lprefix);
                    doDisplay(doc(i_doc).attribute, ntab, dmax, d + 1);
                    fprintf('\n');
                end
                if ~isempty(doc(i_doc).content) && d + 1 <= dmax
                    fprintf('%sContents:\n', lprefix);
                    doDisplay(doc(i_doc).content, ntab, dmax, d + 1);
                    fprintf('\n');
                end
            end
        end
    elseif ischar(doc)
        fprintf('%s%s\n', sprefix, doc);
    elseif ~isempty(doc)
        fprintf('%sUnknown data type\n', sprefix);
    end
end
end

function disptext = fixtext(disptext)
lr = findstr(disptext, '\n');
for i_lr = fliplr(lr)
    disptext = [disptext(1:i_lr) sprefix ...
        disptext((i_lr + 1):end)];
end
end