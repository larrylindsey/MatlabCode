function dump_match_tex(fout, eval_parameters)
% dump_match_tex(fout, eval_parameters)
% Continue the match results' tex file with tables of matching results.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04282008  init code
%

pdf_name_gt = get_storage_file_name([eval_parameters.save_prefix, '.gt.pdf']);
pdf_name_as = get_storage_file_name([eval_parameters.save_prefix, '.as.pdf']);

system(['tiff2pdf ', get_storage_file_name([eval_parameters.save_prefix, '.gt.tif']), ...
  ' > ', pdf_name_gt]);
system(['tiff2pdf ', get_storage_file_name([eval_parameters.save_prefix, '.as.tif']), ...
  ' > ', pdf_name_as]);

fprintf(fout, '\\begin{figure*}[t]\n\\begin{center}\n\\begin{tabular}{c}\n');

fprintf(fout, 'Un/matched manual segmentations\\\\\n');
fprintf(fout, '\\includegraphics[width=10cm,type=pdf,ext=.pdf,read=.pdf]{%s}\\\\\n', ...
  pdf_name_gt(1:end-4));

fprintf(fout, 'Un/matched automatic detections\\\\\n');
fprintf(fout, '\\includegraphics[width=10cm,type=pdf,ext=.pdf,read=.pdf]{%s}\\\\\n', ...
  pdf_name_as(1:end-4));

fprintf(fout, '\\end{tabular}\n\\end{center}\n\\end{figure*}\n');
  
return
end
