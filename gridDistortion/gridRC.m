function rc = gridRC(varargin)

rc = cell(size(varargin));

[rc{:}] = meshgrid(varargin{:});

for ii = 1:numel(rc)
    rcl = rc{ii};
    rc{ii} = rcl(:);
end

rc = cat(2, rc{:});

