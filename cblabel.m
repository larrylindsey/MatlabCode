function chout = cblabel(varargin)
if isnumeric(varargin{1})
    ch = varargin{1};
    intext = varargin{2};
    varargin = varargin(3:end);
else
    ch = colorbar;
    intext = varargin{1};
    varargin = varargin(2:end);
end

chyl = get(ch, 'YLabel');
set(chyl, 'String', intext);

while ~isempty(varargin)
    set(chyl, varargin{1}, varargin{2});
    varargin = varargin(3:end);
end

if nargout > 0
    chout = ch;
end