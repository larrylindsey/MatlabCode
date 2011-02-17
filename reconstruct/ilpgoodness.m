function confidence = ilpgoodness(filename)
% function goodness = ilpgoodness(filename)
%
% Calculates Larry's goodness metric on the output of an ilastik
% label prediction.
%
%   goodness - the goodness of the label prediction.
%
%   filename - path to the ilastik project .ilp or .h5 file.

[prediction label] = readILP(filename);

nl = size(prediction, 1);

ltest = reshape(1:nl, [nl 1]);
ltest = repmat(ltest, size(label));

label_rmvect = ones(ndims(label));
label_rmvect(1) = nl;

labelbw = ltest == repmat(label, label_rmvect);

lct = sum(sum(labelbw, 2), 3);

confidence = sum(sum(labelbw .* prediction, 2), 3) ./ lct;
confidence = prod(confidence);


