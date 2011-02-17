function close_match_tex(fout)
% close_match_tex(fout)
% Close the tex file after all match results have been included.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04292008  init code
%

fprintf(fout, '\\end{document}\n');
fclose(fout);
