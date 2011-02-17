function mitochondria_pipeline_NJS(config, sub_module_id)
% mitochondria_pipeline_NJS(config, sub_module_id)
% Pipeline for mitochondria detection based on Nicholas' work.
%
% Nicholas Sofroniew
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  03122009  init. code.
%

switch(sub_module_id)
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Segment mitochondria from the mother cells.
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 0
    mitochondria_segment_2D(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract features to represent segments of mitochondria to distinguish
    % them from those belonging to rest of the neuropile.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 10
    mitochondria_segment_extract_features(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign ground truth labels to the segments.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 11
    mitochondria_segment_ground_truth_label(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Train a classifier to distinguish betweens segments belonging to
    % mitochondria from others.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 12
    mitochondria_segment_train_classifier(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply trained classifier on images.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 13
    mitochondria_segment_apply_classifier(config);

  case 14
    mitochondria_segment_label_threshold(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute 3D linkage confidences for linking mitochondria segments.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 20
    mitochondria_3D_linkage_confidence(config);
  case 28
    mitochondria_evaluate_3D_linkage(config); % DOES NOT EXIST

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate segment labels
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply Markov Random Field to obtain a joint labelling of segments
    % into mitochondria.
  case 30
    mitochondria_segment_label_MRF(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate detection map from segment labels.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 40
    mitochondria_segment_label_2_detection_map(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post process mitochondria detections. *
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 50
    mitochondria_3D_post_process_segment(config);
  case 55
    mitochondria_post_process_voxel(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate mitochondria detections.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    % Segment-wise
    %%%
  case 70
    mitochondria_evaluate_labeling_segmentwise(config);
    %%%
    % Voxel-wise *
    %%%
  case 71
    mitochondria_evaluate_detection_voxelwise(config);
    % Generate PR curves. *
  case 72
    mitochondria_evaluate_PRcurve(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output routines.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 80
    mitochondria_display_detection(config);
  case 85
    mitochondria_display_3D_volume(config);

end
return
end
