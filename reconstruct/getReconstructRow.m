function [row, csvData] = getReconstructRow(name, csvData)
% function [row csvOut] = getReconstructRow(csvData, csvCell)
% Gets a row from the csv file output from Reconstruct's object menu,
%     corresponding to the given object name.
%
% name - the name of the object in question
% csvData - either a struct as returned by importdata, or the filename of
%     the csv file
%
% row - the row in the csv table corresponding the given object name. If
%       the given name was not found, row is returned empty.
% csvOut - the struct corresponding to the csv file, as returned by
%          importdata.
% 

if ischar(csvData)
    csvData = importdata(csvData);
end

[~, ifound] = strsearch(csvData.textdata(:,1), ['^' name '$']);

if ~isempty(ifound)
    row = csvData.data(ifound,:);
end
