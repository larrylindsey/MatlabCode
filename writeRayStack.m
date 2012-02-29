function writtendata = writeRayStack(data, filename, doDel)

if nargin < 3
    doDel = false;
end

if size(data, 1) > size(data, 3)
    data = permute(data, [3 1 2]);
    disp 'Permuted'
end

if doDel && exist(filename, 'file') > 0
    unix(['rm ' filename]);
end

h5create(filename, '/stack', size(data), 'Datatype', class(data));
h5write(filename, '/stack', data);

writtendata = h5read(filename, '/stack');
