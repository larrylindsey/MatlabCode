clear all
close all

case_ids = {};
case_ids{1} = 311:460;
case_ids{2} = 461:610;
% case_ids{3} = 181:190;
% case_ids{4} = 191:200;
% case_ids{5} = 201:210;
% case_ids{6} = 211:220;
% case_ids{7} = 221:230;
% case_ids{8} = 231:240;
% case_ids{9} = 241:250;
% case_ids{10} = 251:260;
% case_ids{11} = 261:270;
% case_ids{12} = 271:280;
% case_ids{13} = 281:290;
% case_ids{14} = 291:310;

for i = 1:length(case_ids)
  em_reconstruct('shinya_Leginon_311_610_combine_proofread_stack', 750, 'case_ids', case_ids{i});
end
