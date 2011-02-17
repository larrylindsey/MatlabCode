function test_mitochondria()


reconstruction_param.mitochondria.manual.image_dir = '/groups/chklovskii/chklovskiilab/Shiv/Marta_aligned/mitochondria/training/';
reconstruction_param.mitochondria.manual.image_prefix = 'a%04d_crop';
reconstruction_param.mitochondria.manual.image_suffix = '.tif';
reconstruction_param.mitochondria.manual.save_suffix = '_mitochondria.mat';
reconstruction_param.mitochondria.manual.case_ids = [1, 2];

mitochondria_manual_annotation(reconstruction_param);

return

end