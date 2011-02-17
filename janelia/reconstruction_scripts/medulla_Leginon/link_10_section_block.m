clear all
close all

% case_ids = {};
% case_ids{1} = 161:170;
% case_ids{2} = 171:180;
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
% 
% for i = 1:length(case_ids)
%   shinya_Leginon_161_310(750, case_ids{i});
% end

% case_ids = {};
% case_ids{1} = 711:720;
% case_ids{2} = 721:730;
% case_ids{3} = 731:740;
% case_ids{4} = 741:750;
% case_ids{5} = 751:760;
% case_ids{6} = 761:770;
% case_ids{7} = 771:780;
% case_ids{8} = 781:790;
% case_ids{9} = 791:800;
% case_ids{10} = 801:810;
% 
% for i = 1:length(case_ids)
%   em_reconstruct('shinya_Leginon_711_810', 750, 'case_ids', case_ids{i});
% end

% case_ids = {};
% case_ids{1} = 1301:1310;
% case_ids{2} = 1311:1320;
% case_ids{3} = 1321:1330;
% case_ids{4} = 1331:1340;
% case_ids{5} = 1341:1350;
% case_ids{6} = 1351:1360;
% case_ids{1} = 1361:1370;
% 
% for i = 1:length(case_ids)
%   em_reconstruct('shinya_Leginon_1301_1370', 750, 'case_ids', case_ids{i});
% end

% case_ids = {};
% case_ids{1} = 1361:1370;
% case_ids{2} = 1371:1380;
% case_ids{3} = 1381:1390;
% case_ids{4} = 1391:1400;
% case_ids{5} = 1401:1410;
% case_ids{6} = 1411:1420;
% case_ids{7} = 1421:1430;
% case_ids{8} = 1431:1440;
% case_ids{9} = 1441:1450;
% case_ids{10} = 1451:1460;
% 
% for i = 1:length(case_ids)
%   em_reconstruct('shinya_Leginon_1361_1460', 750, 'case_ids', case_ids{i});
% end

case_ids = {};
case_ids{1} = 1406:1410;
case_ids{2} = 1411:1415;
case_ids{3} = 1416:1420;
case_ids{4} = 1421:1425;
case_ids{5} = 1426:1430;
case_ids{6} = 1431:1435;
case_ids{7} = 1436:1440;
case_ids{8} = 1441:1445;
case_ids{9} = 1446:1450;
case_ids{10} = 1451:1455;
case_ids{11} = 1456:1460;

for i = 1:length(case_ids)
  em_reconstruct('shinya_Leginon_1361_1460', 750, 'case_ids', case_ids{i});
end