function seriesFlawSplitter(farg, d)

if nargin < 2
    d = 22.5;
end

imlist = orderImageFiles(farg);

fid = fopen('traklist', 'w');

for ii = 1:numel(imlist)
    fname = imlist{ii};
    
    imsplit = splitImageByFlaws(fname, 2);
    
    if size(imsplit, 3) > 1
        
        splitnames = cell(1, size(imsplit, 3));
        
        idot = find(fname == '.', 1, 'first');
        if isempty(idot)
            idot = numel(fname);
            suffix = '.png';
        else
            suffix = fname(idot:end);
            idot = idot - 1;
        end
        
        prefix = fname(1:idot);
                
        for jj = 1:size(imsplit, 3);
            splitnames{jj} = [prefix sprintf('_split%03d', jj) suffix];
            imwrite(imsplit(:,:,jj), splitnames{jj});
            fprintf(fid, '%s 0 0 %g\n', splitnames{jj}, ii * d);
        end
    else
        fprintf(fid, '%s 0 0 %g\n', fname, ii * d);
    end
    
end

