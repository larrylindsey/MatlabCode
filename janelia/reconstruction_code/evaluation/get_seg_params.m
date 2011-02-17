function seg_params = get_seg_params(segmentation_id)
% seg_params = get_seg_params(segmentation_id)
% Get segmentation parameters for different algorithms
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04282008  init code
%

switch(segmentation_id)
  case 1 % Grayscale, watershed
    seg_params = get_seg_params_grayscale_watershed();
  case 2 % Grayscale, ladder
    seg_params = get_seg_params_grayscale_ladder();
  case 3 % Grayscale hist.eq., ladder
    seg_params = get_seg_params_grayscale_heq_ladder();
  case 4 % pb-1, watershed
    seg_params = get_seg_params_pb_1_watershed();
  case 5 % pb-2, watershed
    seg_params = get_seg_params_pb_2_watershed();
  case 6 % pb-3, watershed
    seg_params = get_seg_params_pb_3_watershed();
  case 7 % pb-4, watershed
    seg_params = get_seg_params_pb_4_watershed();
  case 8 % pb-5, watershed
    seg_params = get_seg_params_pb_5_watershed();
  case 9 % pb-1, ladder
    seg_params = get_seg_params_pb_1_ladder();
  case 10 % pb-2, ladder
    seg_params = get_seg_params_pb_2_ladder();
  case 11 % pb-3, ladder
    seg_params = get_seg_params_pb_3_ladder();
  case 12 % pb-4, ladder
    seg_params = get_seg_params_pb_4_ladder();
  case 13 % pb-5, ladder
    seg_params = get_seg_params_pb_5_ladder();
  case 14 % Grayscale hist.eq., ladder, mitochondria
    seg_params = get_seg_params_grayscale_heq_ladder_mitochondria();
  case 15 % pb-1, ladder, mitochondria
    seg_params = get_seg_params_pb_1_ladder_mitochondria();
  case 16 % pb-2, ladder, mitochondria
    seg_params = get_seg_params_pb_2_ladder_mitochondria();
  case 17 % pb-3, ladder, mitochondria
    seg_params = get_seg_params_pb_3_ladder_mitochondria();
  case 18 % pb-4, ladder, mitochondria
    seg_params = get_seg_params_pb_4_ladder_mitochondria();
  case 19 % pb-5, ladder, mitochondria
    seg_params = get_seg_params_pb_5_ladder_mitochondria();
  case 20 % pb-max over filters, ladder, mitochondria
    seg_params = get_seg_params_pb_max_ladder_mitochondria();
end
return
end
