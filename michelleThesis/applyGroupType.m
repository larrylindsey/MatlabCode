function stat = applyGroupType(stat, animalIdAndGroup)
% stat = applyGroupType(stat, animalIdAndGroup)
%  For an animal statistics struct as generated by michelleStatistics, sets
%  the group type field for each element.
%
%  stat - an animal statistics struct array as generated by
%         michelleStatistics
% animalIDAndGroup - a 2D array in [n 2], where there are n number of
%                    unique animal ids. The first column is taken to be the
%                    animal id, the second the group number.
%
% Typically, the animalID field in stat will be a string, whereas in
% animalIdAndGroup it is numeric. This function handles the disparity
% automagically

n = size(animalIdAndGroup, 1);
m = numel(stat);

idGroupHash = struct;
% use a struct like a hash table
for i_n = 1:n
    animalId = getFKey(animalIdAndGroup(i_n, 1));
    idGroupHash.(animalId) = animalIdAndGroup(i_n, 2);
end

for i_m = 1:m
    animalId = getFKey(stat(i_m).animal);
    stat(i_m).group = idGroupHash.(animalId);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fkey = getFKey(key)
if isnumeric(key) || ischar(key) && strcmp(key, num2str(str2num(key))) == 1 
    fkey = ['n' num2str(key)];
else
    fkey = key;
end
end