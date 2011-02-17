function sandbox()

reconstruction_param.dir = '/groups/chklovskii/chklovskiilab/Shiv/Marta_aligned/';
reconstruction_param.image_prefix = 'a%04d';
reconstruction_param.image_suffix = '.tif';
reconstruction_param.case_ids = 1:20;

roi = get_roi(reconstruction_param)

return

end