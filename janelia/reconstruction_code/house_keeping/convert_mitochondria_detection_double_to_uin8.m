% convert mitochondria detection maps from double to uint8 to save storage
% space.
% value_double is thresholded to lie between -10 and 10.
% value_uint8 = uint*(value_double * 10 + 100).
% value_double = double(value_uint8)/10 - 10.
%

f = ls('-1 *.mitochondria_det_conf.mat');
line_brackets = [0, strfind(f, char(10))];
n_files = length(line_brackets)-1;

for i = 1:n_files
  file_name = f(line_brackets(i)+1:line_brackets(i+1)-1);
  m = load_mitochondria_detection(file_name);
  mitochondria_detection_confidence = m.mitochondria_detection_confidence;
  save_mitochondria_detection(file_name, mitochondria_detection_confidence);
end
