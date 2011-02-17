function mitochondria_NJS_demo(case_ids_train,case_ids_apply)
% mitochondria_NJS_demo(configcase_ids_train,case_ids_apply)
% Runs NJS mitochondria detector from start to finish 
%
% Nicholas Sofroniew,
% Univ. of Cambridge, UK.
% Visitor March 2009, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03162009  init code
%

mitochondria_2007KHLateral_s09_1_91([0,10:12], case_ids_train);
mitochondria_2007KHLateral_s09_1_91([0,10:13,20,30,40,50,55,71,85], case_ids_apply);