function [fd blockhist] = blockHistogramDifference(blockhist)
%

if isstruct(blockhist)
    blockhist = structConvert(blockhist);
end


blockhistROffset = blockhist(2:end, 1:(end-1), :);
blockhistCOffset = blockhist(1:(end - 1), 2:end, :);

blockhistRInt = min(blockhist(1:(end - 1), 1:(end - 1), :), blockhistROffset);
blockhistRInt = sum(blockhistRInt, 3);

blockhistCInt = min(blockhist(1:(end - 1), 1:(end - 1), :), blockhistCOffset);
blockhistCInt = sum(blockhistCInt, 3);

fd = sqrt(blockhistRInt.^2 + blockhistCInt.^2);
end

function blockhist = structConvert(blockhist)

% Guess the struct field name. Generally, expect the struct to have only
% one field, but it still works as long as the histogram is stored in the
% *first* field.
allfields = fieldnames(blockhist); 
sfield = allfields{1};

[blockR, blockC] = size(blockhist);
blockH = numel(blockhist(1).(sfield));


% Slice blockhist into ncpu pieces, where ncpu is the number of cpus
% available to matlab.
ncpu = matlabpool('size');
ubound = round(linspace(numel(blockhist) / ncpu, numel(blockhist), ncpu));
if ncpu > 1
    lbound = [1 ubound(1:(end - 1)) + 1];
else
    lbound = 1;
end

slice = cell(1, ncpu);
catblockhist_c = slice;

for ii = 1:ncpu
    slice{ii} = blockhist(lbound(ii):ubound(ii));
end

% reshape the block struct into a matlab array
parfor ii = 1:ncpu
    catblockhist_c{ii} = cat(1, slice{ii}.(sfield));
end

blockhist = cat(1, catblockhist_c{:});

blockhist = reshape(blockhist, [blockR, blockC, blockH]);
end
