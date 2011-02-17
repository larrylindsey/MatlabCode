
function tmp()

iids = imgList('test');
ignore = mkdir2('pb');
ignore = mkdir2('pb/bgtg');
for iid = iids,
  im = rgb2gray(imgRead(iid));
  fprintf(2,'Computing Pb for image %d using BG/TG...\n',iid);
  pb = pbBGTG(im);
  imwrite(pb,sprintf('pb/bgtg/%d.bmp',iid),'bmp');
end

