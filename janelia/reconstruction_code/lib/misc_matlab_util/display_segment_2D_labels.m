function d = display_segment_2D_labels(label_map, mapping, figure_id)

if(nargin<2)
  mapping = [];
end
if(nargin<3)
  figure_id = [];
end

if(isempty(figure_id))
  figure;
else
  figure(figure_id);
end

if(~isempty(mapping))
  label_map = apply_mapping(label_map, mapping);
end

label_map = relabel_to_remove_nonexistent_labels(label_map);

mapping_rand = [0 randperm(max(label_map(:)))];
label_map = mapping_rand(label_map+1);

cm = [ ...
         0         0    0
         0         0    0.5625
         0         0    0.6250
         0         0    0.6875
         0         0    0.7500
         0         0    0.8125
         0         0    0.8750
         0         0    0.9375
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0625    1.0000    0.9375
    0.1250    1.0000    0.8750
    0.1875    1.0000    0.8125
    0.2500    1.0000    0.7500
    0.3125    1.0000    0.6875
    0.3750    1.0000    0.6250
    0.4375    1.0000    0.5625
    0.5000    1.0000    0.5000
    0.5625    1.0000    0.4375
    0.6250    1.0000    0.3750
    0.6875    1.0000    0.3125
    0.7500    1.0000    0.2500
    0.8125    1.0000    0.1875
    0.8750    1.0000    0.1250
    0.9375    1.0000    0.0625
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0
    0.9375         0         0
    0.8750         0         0
    0.8125         0         0
    0.7500         0         0
    0.6875         0         0
    0.6250         0         0
    0.5625         0         0
    0.5000         0         0
];

label_map_mod = label_map;
label_map_mod(label_map_mod>0) = mod(label_map_mod(label_map_mod>0), 64)+1;

label_map_mod2 = label_map;
label_map_mod2(label_map_mod2>0) = ...
  mod(label_map_mod2(label_map_mod2>0), 64)+1;

d = zeros([size(label_map), 3]);
for c = 1:3
  cmc = cm(:,c);
  d(:,:,c) = (cmc(label_map_mod+1) + cmc(label_map_mod2+1))/2;
end
d = uint8(255*d);
imshow(d);

return
end
