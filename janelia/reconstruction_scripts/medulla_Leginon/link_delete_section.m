clear all
close all

% case_ids = [768, 779, 788, 793, 799];
% 
% for i = 1:length(case_ids)
%   em_reconstruct('shinya_Leginon_711_810', 750, 'case_ids', case_ids(i)-1:case_ids(i)+1);
% end

% case_ids = [920, 940, 960, 980, 1000];
% 
% for i = 1:length(case_ids)
%   em_reconstruct('shinya_Leginon_911_1010', 750, 'case_ids', case_ids(i):case_ids(i)+1);
% end

% case_ids = [1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090, 1100];
% 
% for i = 1:length(case_ids)
%   em_reconstruct('shinya_Leginon_1011_1110', 750, 'case_ids', case_ids(i):case_ids(i)+1);
% end

case_ids = [1120, 1130, 1140, 1150];

for i = 1:length(case_ids)
  em_reconstruct('shinya_Leginon_1111_1154', 750, 'case_ids', case_ids(i):case_ids(i)+1);
end