function [imstack ids] = read_sop_stack(sz)

dd = dir('*.png');
imstack = zeros([sz numel(dd)]);
ids = zeros(numel(dd), 1);
for ii = 1:numel(dd)
    [imstack(:, :, ii) ids(ii)] = read_sop_image(dd(ii).name, sz);
end

end
