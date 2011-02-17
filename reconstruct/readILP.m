function [prediction label data] = readILP(filename)
% function [prediction label data] = readILP(filename)
%
% Makes reading an ilastik project file a little more convenient.
%
%   filename - path to an ilastik .ilp or .h5 file.
%
%   prediction - label prediction, computed by ilastik
%   label - annotations used for training by ilastik
%   data - the image data used in the project file
%

% This function is pretty dumb - we assume that all of our data is stored
% in dataItem00.  This may or may not work for a 3D project, and it
% definitely won't work if you do an edit->add stack in ilastik.

prediction = hdf5read(filename, '/DataSets/dataItem00/prediction');

if nargout > 1
    label = hdf5read(filename, '/DataSets/dataItem00/labels');
    if nargout > 2
        data = hdf5read(filename, '/DataSets/dataItem00/data');
    end
end


end