function outstr = strsearch(strs, pat)

in = regexp(strs, pat);
outstr = {};

for ii = 1:numel(in)
    if ~isempty(in{ii})
        outstr = {outstr{:}, strs{ii}};
    end
end
