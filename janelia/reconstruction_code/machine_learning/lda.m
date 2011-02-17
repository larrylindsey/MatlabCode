function filters = lda(patches, labels, shape)
    if ~exist('shape', 'var'),
        shape = size(patches, 2);
    end
    [~, ~, stats] = manova1(patches, labels);
    sig_cols = 1:(numel(unique(labels))-1);
    filters = stats.eigenvec(:,sig_cols);
    if size(filters, 2) > 1,
        shape = [shape size(filters, 2)];
    end
    filters = reshape(filters, shape);
end