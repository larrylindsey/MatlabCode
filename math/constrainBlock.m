function trconstrblock = constrainBlock(trblock)
trblock = trblock(:);
trconstrblock = constrainTransform(trblock(1));
trconstrblock = repmat(trconstrblock, [numel(trblock), 3]);

for ii = 1:numel(trblock)
    trconstrblock(ii, 1) = constrainTransform(trblock(ii));
    trconstrblock(ii, 2) = constrainTransform2(trblock(ii));
    trconstrblock(ii, 3) = constrainTransformIxy(trblock(ii));
end
