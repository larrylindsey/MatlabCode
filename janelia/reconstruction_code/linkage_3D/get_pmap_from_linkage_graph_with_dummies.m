function pmap = get_pmap_from_linkage_graph_with_dummies(links_3D, link_threshold)
% pmap = get_pmap_from_linkage_graph(links_3D, link_threshold)
% Compute pmap from the linkage graph by thresholding link-edge weights at
% the link_threshold. Uses a graph parsing algorithm.
%
% Linkage graph has nodes corresponding to segments and links corresponding
% to 3D linkage likelihoods. Connected components in the thresholded graph
% correspond to bodies in the reconstruction. They are obtained using a
% simple graph parsing algorithm.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~03312008  init code
%

n_seg = max(max(links_3D{1}(:,1:2)));
min_dummy_seg_id = 0;
for i = 1:length(links_3D)
  if(isempty(links_3D{i}))
    continue;
  end
  n_seg = max(n_seg, max(max(links_3D{i}(:,1:2))));
  min_dummy_seg_id = min(min_dummy_seg_id, min(min(links_3D{i}(:, 1:2))));
end

pmap = 1:n_seg;

processed_flag = zeros(n_seg, 1);

processed_flag_dummy = zeros(-min_dummy_seg_id,1);

fprintf('building pmap from linkage graph\n');
for i = 1:length(links_3D)
  fprintf('%d ', i);
  for j = 1:size(links_3D{i}, 1)
    if(links_3D{i}(j,1)<=0 || links_3D{i}(j,2)==0)
      continue;
    end
    
    adam = links_3D{i}(j,1);
    if(processed_flag(adam)==1)
      continue;
    end;
    processed_flag(adam)=1;
    
    link_seg = links_3D{i}(links_3D{i}(:,1)==adam & ...
      links_3D{i}(:,3)>link_threshold, 2);
    
    link_list = [ones(length(link_seg), 1)*i+1, link_seg];
    tail = size(link_list,1);
    while(tail>0)
      z = link_list(tail,1);
      seg_id = link_list(tail,2);
      tail = tail-1;
      
      if(seg_id>0) % valid segment id
        if(processed_flag(seg_id)==1)
          continue;
        end;
        processed_flag(seg_id)=1;
        
        pmap(pmap==seg_id) = adam;
      else
        if(processed_flag_dummy(-seg_id)==1)
          continue;
        end;
        processed_flag_dummy(-seg_id)=1;
      end
      
      if(z>1)
        if(~isempty(links_3D{z-1}))
          link_seg = links_3D{z-1}(links_3D{z-1}(:,2)==seg_id  & ...
            links_3D{z-1}(:,3)>link_threshold, 1);
          link_seg_valid = link_seg(link_seg>0);
          link_seg_dummy = link_seg(link_seg<0);
          link_seg_valid = link_seg_valid(processed_flag(link_seg_valid)==0);
          link_seg_dummy = ...
            link_seg_dummy(processed_flag_dummy(-link_seg_dummy)==0);
          link_seg = [link_seg_valid; link_seg_dummy];
          link_list = [link_list; ones(length(link_seg), 1)*z-1, link_seg]; %#ok<AGROW>
          tail = size(link_list,1);
        end
      end;
      if(z<=length(links_3D))
        if(~isempty(links_3D{z}))
          link_seg = links_3D{z}(links_3D{z}(:,1)==seg_id  & ...
            links_3D{z}(:,3)>link_threshold, 2);
          link_seg_valid = link_seg(link_seg>0);
          link_seg_dummy = link_seg(link_seg<0);
          link_seg_valid = link_seg_valid(processed_flag(link_seg_valid)==0);
          link_seg_dummy = ...
            link_seg_dummy(processed_flag_dummy(-link_seg_dummy)==0);
          link_seg = [link_seg_valid; link_seg_dummy];
          link_list = [link_list; ones(length(link_seg), 1)*z+1, link_seg]; %#ok<AGROW>
          tail = size(link_list,1);
        end
      end;
    end;
  end;
end;
fprintf('done\n');
pmap = uint32([0, pmap]);
return;
