function varargout = inputStruct(varargin)

arg1 = varargin{1};

if ischar(arg1)
    outStr = struct;
    for i_arg = 1:nargin
        outStr.(varargin{i_arg}) = [];
    end
    varargout{1} = outStr;
else
    inStr = varargin{1};
    for i_arg = 2:nargin
        varargout{i_arg - 1} = inStr.(varargin{i_arg});
    end
end

