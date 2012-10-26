function testExtraction(im, K, n, rootName)

pk = .5.^((1:n) - 1);

kk = zeros(1, 2*n+1);
kk(1:n) = (-K)*pk;
kk(end:-1:(n+2)) = K*pk;
kk(n+1) = 0;

for i_k = 1:numel(kk)
    fprintf('Working on %d\n', i_k);
    fprintf('\n');
    outname = sprintf('%s%d.mat', rootName, i_k);
    fprintf('Eventually saving output to %s\n\n', outname);
    imGD = quickImDistort(im, kk(i_k), .5);
    imshow(imGD);
    drawnow;
    [tr, db, e1] = tryExtract(imGD);
    imTr = imtransform(imGD, tr);
    imshow(imTr);
    [tr(2), db(2), e2] = tryExtract(imTr);
    fprintf('\n\n');
    fprintf('Saving output now...\n');
    save(outname, 'tr', 'db', 'K', 'e1', 'e2');
    clear outstr;
    fprintf('Done.\n\n');
end

end


function [tr db e] = tryExtract(im)
badcount = 0;
itDidntWork = true;
while itDidntWork && badcount < 2
    try
        [tr, db] = extractDistortionTransform(im);
        itDidntWork = false;
        e = [];
    catch
        e = cat(1, e, lasterror);
        badcount = badcount + 1;
        fprintf('Detected Error, count %d\n', badcount);
    end
end

end