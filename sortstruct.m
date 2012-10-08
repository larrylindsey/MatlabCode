function structout = sortstruct(structin, field)

if ischar(structin(1).(field))
    sortv = {structin.(field)};
else
    sortv = [structin.(field)];
end

[~, isort] = sort(sortv);
structout = structin(isort);
