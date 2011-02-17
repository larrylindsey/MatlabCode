load ../medulla.GA5min_HPF_FS_0509.5000x.000_002_2/linkage_3D/link_dump_k_fold.boost.intensity_pair_hist_v2c.05202008duk.mat

%%
[junk, sort_id] = sort([link_dump.to_be_linked(:).link_conf]);
link_dump.to_be_linked = link_dump.to_be_linked(sort_id);

[junk, sort_id] = sort([link_dump.not_to_be_linked(:).link_conf], 2, 'descend');
link_dump.not_to_be_linked = link_dump.not_to_be_linked(sort_id);

%%
ind = find([link_dump.to_be_linked(:).z]==1);
for i = ind
  figure(100);
  d = repmat(link_dump.to_be_linked(i).image_1/255, [1 1 3]);
  d(:,:,3) = link_dump.to_be_linked(i).segment_1/255;
  subplot(1,2,1); imshow(d, []);
  d = repmat(link_dump.to_be_linked(i).image_2/255, [1 1 3]);
  d(:,:,3) = link_dump.to_be_linked(i).segment_2/255;
  subplot(1,2,2); imshow(d, []);
  title([num2str(link_dump.to_be_linked(i).link_conf), ', ', ...
    num2str(size(d,1)), ', ', num2str(size(d,2))]);
  pause;
end

%%
ind = find([link_dump.not_to_be_linked(:).z]==1);
for i = ind
  figure(101);
  d = repmat(link_dump.not_to_be_linked(i).image_1/255, [1 1 3]);
  d(:,:,3) = link_dump.not_to_be_linked(i).segment_1/255;
  subplot(1,2,1); imshow(d, []);
  d = repmat(link_dump.not_to_be_linked(i).image_2/255, [1 1 3]);
  d(:,:,3) = link_dump.not_to_be_linked(i).segment_2/255;
  subplot(1,2,2); imshow(d, []);
  title([num2str(i), ':', num2str(link_dump.not_to_be_linked(i).link_conf), ', ', ...
    num2str(size(d,1)), ', ', num2str(size(d,2))]);
  pause;
end

%%

ind = find([link_dump.not_to_be_linked(:).link_conf]>0);
for i = ind
  figure(101);
  d = repmat(link_dump.not_to_be_linked(i).image_1/255, [1 1 3]);
  d(:,:,3) = link_dump.not_to_be_linked(i).segment_1/255;
  subplot(1,3,1); imshow(d, []);
  d = repmat(link_dump.not_to_be_linked(i).image_2/255, [1 1 3]);
  d(:,:,3) = link_dump.not_to_be_linked(i).segment_2/255;
  subplot(1,3,2); imshow(d, []);
  d = repmat(1-link_dump.not_to_be_linked(i).image_1/255, [1 1 3]);
  d(:,:,3) = 1-link_dump.not_to_be_linked(i).image_2/255;
  subplot(1,3,3); imshow(d, []);
  title([num2str(link_dump.not_to_be_linked(i).link_conf), ', ', ...
    num2str(size(d,1)), ', ', num2str(size(d,2))]);
  pause;
end

%%

ind = find([link_dump.to_be_linked(:).link_conf]<0);
for i = ind
  figure(100);
  d = repmat(link_dump.to_be_linked(i).image_1/255, [1 1 3]);
  d(:,:,3) = link_dump.to_be_linked(i).segment_1/255;
  subplot(1,3,1); imshow(d, []);
  d = repmat(link_dump.to_be_linked(i).image_2/255, [1 1 3]);
  d(:,:,3) = link_dump.to_be_linked(i).segment_2/255;
  subplot(1,3,2); imshow(d, []);
  d = repmat(1-link_dump.to_be_linked(i).image_1/255, [1 1 3]);
  d(:,:,3) = 1-link_dump.to_be_linked(i).image_2/255;
  subplot(1,3,3); imshow(d, []);
  title([num2str(link_dump.to_be_linked(i).link_conf), ', ', ...
    num2str(size(d,1)), ', ', num2str(size(d,2))]);
  pause;
end
