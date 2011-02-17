function submit_to_cluster_copy_global_patchwise_transforms(case_ids, jobs_dir, ...
  bin_dir, submit_script_name)

config = get_basic_config();

config.stack.dir = '/groups/chklovskii/medulla/';
config.stack.name = 'medulla.HPF.Leginon.3500x.zhiyuan.fall2008';   
config.stack.image_structure = 'crop4_global_alignment_0161_1460.kludge.unreal.xml';
config.stack.image_structure_is_kludged = true;
fprintf('case_ids: %d\n', case_ids);
config.stack.case_ids = case_ids;
config.region.case_ids = 161:1460;
config.stack.roi = [];

region_structure = get_region_structure(config);

fout = fopen([jobs_dir, submit_script_name], 'wt');

z = region_structure.planes(1).z;
output_dir = [jobs_dir, num2str(z), '/'];
check_for_dir(output_dir);
save([output_dir, 'region_structure.mat'], '-STRUCT', 'region_structure');

fprintf(fout, 'cd %d\n', z);
fprintf(fout, ['qsub -N c%d -l excl=true -j y -o log -b y -cwd ', ...
  '-V "%scopy_global_patchwise_transforms_from_simple_txt_to_mat_files ', ...
  'region_structure.mat"\n'], z, bin_dir);
fprintf(fout, 'cd ..\n');

fclose(fout);

return
end
