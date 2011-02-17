function pmap = get_pmap_from_linkage_graph(links_3D, link_threshold)
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
for i = 1:length(links_3D)
  n_seg = max(n_seg, max(max(links_3D{i}(:,1:2))));
end

pmap = 1:n_seg;

processed_flag = zeros(n_seg, 1);

fprintf('building pmap from linkage graph\n');
for i = 1:length(links_3D)
  fprintf('%d ', i);
  for j = 1:size(links_3D{i}, 1)
    if(min(links_3D{i}(j,1:2))<=0)
      continue;
    end
    
    adam = links_3D{i}(j,1);
    if(processed_flag(adam)==1)
      continue;
    end;
    processed_flag(adam)=1;
    
    link_seg = links_3D{i}(links_3D{i}(:,1)==adam & links_3D{i}(:,3)>link_threshold, 2);
    
    link_list = [ones(length(link_seg), 1)*i+1, link_seg];
    tail = size(link_list,1);
    while(tail>0)
      z = link_list(tail,1);
      seg_id = link_list(tail,2);
      tail = tail-1;
      
      if(processed_flag(seg_id)==1)
        continue;
      end;
      processed_flag(seg_id)=1;
      
      pmap(pmap==seg_id) = adam;
      
      if(z>1)
        link_seg = links_3D{z-1}(links_3D{z-1}(:,2)==seg_id  & links_3D{z-1}(:,3)>link_threshold, 1);
        link_seg = link_seg(processed_flag(link_seg)==0);
        link_list(tail+1:tail+length(link_seg), :) = [ones(length(link_seg), 1)*z-1, link_seg];
        tail = size(link_list,1);
      end;
      if(z<=length(links_3D))
        link_seg = links_3D{z}(links_3D{z}(:,1)==seg_id  & links_3D{z}(:,3)>link_threshold, 2);
        link_seg = link_seg(processed_flag(link_seg)==0);
        link_list(tail+1:tail+length(link_seg), :) = [ones(length(link_seg), 1)*z+1, link_seg];
        tail = size(link_list,1);
      end;
      
    end;
    
  end;
end;
fprintf('done\n');
pmap = uint32([0, pmap]);
return;