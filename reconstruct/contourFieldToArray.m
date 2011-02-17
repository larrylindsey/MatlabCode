function output = contourFieldToArray(strarray, field)

output = cell(0);

for i_str = 1:numel(strarray)
    fieldval = strarray(i_str).(field);
    output = {output{:} fieldval};
end