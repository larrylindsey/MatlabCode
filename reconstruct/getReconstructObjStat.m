function [stat, names, csvData] = getReconstructObjStat(pattern, csvData, exact)
% function [row names csvOut] = getReconstructRow(pattern, csvData, exact)
% Gets a row from the csv file output from Reconstruct's object menu,
%     corresponding to the given object name.
%
% pattern - the pattern matching the object(s) in question
% csvData - either a struct as returned by importdata, or the filename of
%     the csv file
% exact - true if the object's name should match pattern exactly, false to
%         match as a regular expression
%
% stat - a struct containing the measurements taken from the row in the csv
%        file corresponding to the given object name
% names - a cell array containing the matched names from the csvData struct
% csvOut - the struct corresponding to the csv file, as returned by
%          importdata.
% 

if ischar(csvData)
    csvData = importdata(csvData);
end

if nargin < 3
    exact = false;
end

if exact
    [~, ifound] = strsearch(csvData.textdata(2:end,1), ['^' pattern '$']);
    names = {pattern};
else
    [names, ifound] = strsearch(csvData.textdata(2:end,1), pattern);
end

if ~isempty(ifound)
    row = csvData.data(ifound,:);
end

stat = repmat(struct(), [numel(ifound) 1]);

for jj = 1:numel(ifound)
    stat(jj).Object = names{jj};
    for ii = 1:size(row,2)
        tf = csvData.textdata{1,ii + 1};
        tf(tf == ' ') = '';
        stat(jj).(tf) = row(jj,ii);        
    end    
end

stat = sortstruct(stat, 'Object');
