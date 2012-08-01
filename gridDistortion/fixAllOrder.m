function eanalysis = fixAllOrder(eanalysis, dirs)

f = figure;

for ii = 1:numel(dirs)
    imfile = ls([dirs{ii} '/*.0.jpg']);
    imfile(end) = ''; 
    im = imread(imfile);
    im = im2single(im(:,:,1));
    im = imresize(im, .25);
    
    jj = 1;
    trblock = eanalysis(ii).trblock;
    iorder = eanalysis(ii).allorder;
    iok = [];
    while jj < numel(iorder) && numel(iok) < 5
        imtr = applyTransformImage(im, trblock(iorder(jj)));
        imshow(imtr);
        figure(f);
        p = ginput(1);
        if ~isempty(p)
            iok = [iok jj];
            disp('Accepted!');
        else
            disp('Rejected!');
        end
        jj = jj + 1;
    end
    eanalysis(ii).iok = iok;    
end
