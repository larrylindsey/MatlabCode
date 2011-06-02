function xmlCell = loadXML(file)

fid = fopen(file, 'r');
if fid < 1
    error('Failed to open file %s', file);
end

currLine = fgetl(fid);


il = 1;

while ischar(currLine)
    tok = strtok(currLine);
    if ~isempty(tok)
        xmlCell{il} = currLine;
        il = il + 1;
    end
    currLine = fgetl(fid);
end

if ~ischar(currLine)
    if il == 1
        xmlCell = {};
    end
end



fclose(fid);