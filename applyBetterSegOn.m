function applyBetterSegOn(liststr)

imlist = dir(liststr);
imnames = {imlist.name};

disp(numel(imnames));

for ii = 1:numel(imnames)
    filename = imnames{ii};
    segname = ['class_' filename];
    fprintf('%d, reading files\n', ii);
    im = rgb2gray(imread(filename));
    seg = 255 - imread(segname);
    diffsize = (size(seg) - size(im)) + 1;
    diffsize = floor(diffsize / 2);
    seg = seg(diffsize(1):(end - diffsize(1)), diffsize(2):(end - diffsize(2)));
    seg = seg(1:size(im,1), 1:size(im,2));
    fprintf('%d, fixing up segment\n', ii);
    betterseg = better_segment(seg, im) / 65536;
    fprintf('%d, writing out segmented image\n', ii)
    imwrite(betterseg, ['seg_' filename '.png'], 'PNG', 'BitDepth', 16);
    fprintf('%d done\n', ii);
end