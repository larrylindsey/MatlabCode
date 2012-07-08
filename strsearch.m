function [outstr ifound] = strsearch(strs, pat)

in = regexp(strs, pat);
outstr = {};
ifound = [];

for ii = 1:numel(in)
    if ~isempty(in{ii})
        outstr = [outstr, strs(ii)]; %#ok<*AGROW>
        ifound = [ifound ii];
    end
end
