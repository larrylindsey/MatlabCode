function varargout = errorEmail(func, varargin)

if ischar(func)
    func = str2func(func);
end

varargout = cell(1, nargout);
try
[varargout{:}] = func(varargin{:});
catch e
    [junk hname] = unix('hostname');
    cc = clock;
    tstr = sprintf('%02d-%02d-%4d %02d:%02d:%02d', cc(3), cc(2), cc(1),...
        cc(4), cc(5), round(cc(6)));
    str = sprintf('Caught error on %s at %s\n%s\nIn %s in %s\nAt line %d', ...
        hname, tstr, e.message, e.stack(1).name, e.stack(1).file, e.stack(1).line);
    sendmail('larry.f.lindsey@gmail.com','Caught Error', str);
    %rethrow(e);
end
end