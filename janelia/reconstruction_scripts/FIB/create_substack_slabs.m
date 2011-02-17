% create substacks, e.g., of 30 sections for 3D segmentation
for z = 1:30:900
    fprintf('z: %d\n', z);
    em_reconstruct('recon3D_lamina_FIB_5x5x5nm_esb_denoised', 1001, 'case_ids', z:z+29);
end

start_indexes = 1:30:900;
end_indexes = start_indexes+29;
i = [start_indexes; end_indexes];
fprintf('\n');
fprintf('%d %d ', i(:));
fprintf('\n');
