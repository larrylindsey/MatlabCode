function save_mitochondria_detection(file_name, mitochondria_detection_confidence)
% Wrapper for saving mitochondria detections

mitochondria_detection_confidence = mitochondria_detection_confidence*10 + 100;

mitochondria_detection_confidence(mitochondria_detection_confidence>200) = 200;
mitochondria_detection_confidence(mitochondria_detection_confidence<0) = 0;


mitochondria_detection_confidence = uint8(mitochondria_detection_confidence); %#ok<NASGU>

save2(file_name, 'mitochondria_detection_confidence');

return
end
