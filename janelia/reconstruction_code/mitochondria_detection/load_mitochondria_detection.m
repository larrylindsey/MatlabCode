function mitochondria_detect = load_mitochondria_detection(file_name)
% Wrapper for loading mitochondria detections

mitochondria_detect = load2(file_name);
if(strcmp(class(mitochondria_detect.mitochondria_detection_confidence), 'uint8')==1)
  mitochondria_detect.mitochondria_detection_confidence = ...
    double(mitochondria_detect.mitochondria_detection_confidence);
  mitochondria_detect.mitochondria_detection_confidence = ...
    mitochondria_detect.mitochondria_detection_confidence/10;
  mitochondria_detect.mitochondria_detection_confidence = ...
    mitochondria_detect.mitochondria_detection_confidence - 10;
end

return
end
