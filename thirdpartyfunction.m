function varargout = thirdpartyfunction(func, varargin)

codedir = '/home/larry/code/3rd-party-matlab';

if ischar(func)
    func = str2func(func);
end

funcid = functions(func);
if isempty(funcid.file)
    cmd = sprintf(['find %s/ -name %s.* '...
        '| sed -e "%s" | gawk ''{ print $5 }'''], ...
        codedir, funcid.function, 's/\// /g');
    [status fdir] = unix(cmd);
    
    if status ~= 0
        error('Encountered unix error: %s', fdir);
    end
    
    fdir = textscan(fdir, '%s');
    adir = fdir{1};
    
    for ii = 1:numel(adir)
        addpath(genpath([codedir '/' adir{ii}]));
    end
end

varargout = cell(1, nargout);

[varargout{:}] = func(varargin{:});
end