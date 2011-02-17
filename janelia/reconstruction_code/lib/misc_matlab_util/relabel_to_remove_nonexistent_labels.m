function [s_relabelled, relabel_map] = relabel_to_remove_nonexistent_labels(s)
% [s_relabelled, relabel_map] = relabel_to_remove_nonexistent_labels(s)
% Relabels a matrix of labels, 's', to remove non-existent labels.
% Labels are assumed to be non-negative.
% Outputs:
%   s_relabelled    relabelled matrix
%   relabel_map     mapping such that s_relabelled = relabel_map(1+s);
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

label_ids = [0; nonzeros(unique(s(:)))];
relabel_map = [];
relabel_map(1+label_ids) = 0:(length(label_ids)-1);
s_relabelled = relabel_map(1+s);

return
end
