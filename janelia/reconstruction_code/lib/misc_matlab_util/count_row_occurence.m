function [rows_unique, count] = count_row_occurence(s, is_sorted)
% [rows_unique, count] = count_row_occurence(s, is_sorted)
% Count the number of times each row occurs in s. Optionally, set is_sorted
% to true if s is pre-sorted for efficiency.
%
% Output:
%   rows_unique     unique rows in s
%   count           their counts: count(i) is the number of times
%                     rows_unique(i,:) occurs in s.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(nargin<2)
  is_sorted = false;
end

if(~is_sorted)
  s_sorted = sortrows(s);
else
  s_sorted = s;
end

[rows_unique, id_first] = unique(s_sorted, 'rows', 'first');
[junk, id_last] = unique(s_sorted, 'rows', 'last');
count = id_last - id_first + 1;

return;
end
