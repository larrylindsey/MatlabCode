function transformFiles(tr, infiles)
for i_f = 1:numel(infiles)
    file = infiles{i_f};
    outfile = [file(1:(end-4)) 'trans' file((end-3):end)];
    im = imread(file);
    imTrans = imtransform(im, tr);
    fprintf('%s -> %s.  %g of %g\n', file, outfile, i_f, numel(infiles));    
    imwrite(imTrans, outfile);
end