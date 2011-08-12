function list = orderImageFiles(darg, p)
% function list = orderImageFiles(darg, [p])
%
% Sorts a list of files by correct numerical order.
%  darg - may be an argument to dir, such as '*.jpg', or a cell array of 
%         file names.
%
%  p - if darg is not a cell array, p should be the directory that the 
%      files reside in. If it is omitted, p takes the value returned by
%      pwd.
%
%  list - A sorted cell array of the file names in question.


if ~iscell(darg)
    if nargin < 2
        p = pwd;
    end
    
    d = dir([p '/' darg]);
    
    prelist = {d.name};
    
else
    prelist = darg;
end

prefix = prelist{1};
suffix = prelist{1};
suffix = suffix(end:-1:1);
ns = numel(suffix);
np = ns;

for ii = 2:numel(prelist)
    currname = prelist{ii};
    ns = min(ns, stringsAlignTo(suffix, currname(end:-1:1)));
    np = min(np, stringsAlignTo(prefix, currname));
end

sortlist = zeros(size(prelist));

for ii = 1:numel(prelist)
    sortstr = prelist{ii}((1 + np):(end - ns));
    sortlist(ii) = getNumericValue(sortstr);
end

[junk order] = sort(sortlist, 'ascend'); %#ok

list = prelist(order);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = stringsAlignTo(str1, str2)
l = min(numel(str1), numel(str2));

n = find(str1(1:l) ~= str2(1:l), 1, 'first') - 1;

if isempty(n)
    n = 0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = getNumericValue(str)

strN = uint8(str);
sel = and(strN >= 48, strN <= 57);

s = find(sel, 1, 'first');

if isempty(s)
    v = -1;
else
    l = find(not(sel(s:end)), 1, 'first') - 1;
    if isempty(l)
        l = numel(sel(s:end));
    end
    v = str2double(str((1:l) - 1 + s));
end

end

