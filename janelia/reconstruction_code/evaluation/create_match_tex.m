function create_match_tex(fout, title_text)
% create_match_tex(fout)
% Put the bipartite match images into a tex for compilation to pdf using
% pdflatex. create_match_pdf() starts the preamble of the tex file. To be
% continued with tables of matching results.
%
fprintf(fout, '\\documentclass{article}\n\\usepackage{graphicx}\n');
fprintf(fout, '\\begin{document}\n\\title{Results of Evaluating Segmentation ');
fprintf(fout, 'using CSA++ Algorithm}\n\\author{}\n\\maketitle\n');
title_text = strrep(title_text, '_', ' ');
fprintf(fout, 'Approach: %s\n', title_text);
return
end
