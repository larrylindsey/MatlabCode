function testExtractionNoisy(im, K, n, rootName, s)

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
    imGD = addNoise(imGD, s);
    imshow(imGD);
    drawnow;
    [tr1, db1, e1] = tryExtract(imGD);
    imTr = imtransform(imGD, tr1, 'FillValues', .5);
    imshow(imTr);
    [tr2, db2, e2] = tryExtract(imTr);
    fprintf('\n\n');
    fprintf('Saving output now...\n');
    tr = {tr1, tr2};
    db = {db1, db2};
    e = {e1, e2};
    save(outname, 'tr', 'db', 'K', 'e');
    clear outstr;
    fprintf('Done.\n\n');
end

end


function [tr db e] = tryExtract(im)
badcount = 0;
itDidntWork = true;
e = [];
tr = struct;
db = struct;

while itDidntWork && badcount < 2
    try       
        [tr, db] = extractDistortionTransform(im);
        itDidntWork = false;
    catch
        e = cat(1, e, lasterror);
        badcount = badcount + 1;
        fprintf('Detected Error, count %d\n', badcount);
    end
end

end

function im = addNoise(im, s)
noise = randn(size(im)) * s;
im = im + noise;
im(logical(im < 0)) = 0;
im(logical(im > 1)) = 1;
end