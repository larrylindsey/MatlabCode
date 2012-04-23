function varargout = thirdpartyfunction(func, varargin)
% <[argout1 argout2 ... argoutM]> =
%                      thirdpartyfunction(func, <argin1 argin2 ... arginN>)
%  
%  Finds and executes the given function, when it is not part of the path,
%   but exists in a pre-defined directory.
%
% func - the function to execute. May be either a string or a function
%        handle
%
%
% This function requires the 'awk' executable to be present and in
% Matlab's execution path.


codedir = '/home/larry/code/3rd-party-matlab';

dohelp = false;

if strcmp(func, 'help')
    dohelp = true;
    func = varargin{1};
end

if ischar(func)    
    func = str2func(func);
end

funcid = functions(func);
if isempty(funcid.file)
    cmd = sprintf(['find %s/ -name %s.* '...
        '| sed -e "%s" | awk ''{ print $5 }'''], ...
        codedir, funcid.function, 's/\// /g');
    [status fdir] = unix(cmd);
    
    if status ~= 0
        error('Encountered unix error: %s', fdir);
    end
    
    if isempty(fdir)
        error('Function %s not found in %s', funcid.function, codedir);
    end
    
    fdir = textscan(fdir, '%s');
    adir = fdir{1};
    
    for ii = 1:numel(adir)
        addpath(genpath([codedir '/' adir{ii}]));
    end
end

varargout = cell(1, nargout);
if dohelp
    help(funcid.function);
else    
    [varargout{:}] = func(varargin{:});
end
end
