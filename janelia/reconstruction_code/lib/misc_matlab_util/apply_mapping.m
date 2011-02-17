function labels_mapped = apply_mapping(labels, mapping)
% m = [];
% m(1+mapping(:,1)) = mapping(:,2);
% labels_mapped = m(1+labels);

m = zeros([1+max(labels(:)),1]);
m(1+mapping(:,1)) = mapping(:,2);
labels_mapped = m(1+labels);

return
end
