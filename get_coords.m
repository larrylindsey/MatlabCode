function coords = get_coords(basis, el, xx, yy)


sz_el = size(el);
sz_b = size(basis);

if nargin < 3
    [xx yy] = meshgrid(linspace(-1,1,sz_el(1)), linspace(-1,1,sz_el(2)));
end

if sz_el(1:2) ~= sz_b(1:2)
    error('element and basis must have the same support');
end

m_plicity = prod(sz_el(3:end));
coords = zeros(m_plicity, sz_b(3));

el = reshape(el, [sz_el(1:2), m_plicity]);

for i_b = 1:sz_b(3)
    bDot = repmat(basis(:,:,i_b),[1,1,m_plicity]) .* el;
    for i_m = 1:m_plicity
        coords(i_m, i_b) = trapint2(xx,yy,bDot(:,:,i_m));
    end
end

if numel(sz_el(3:end)) == 0
    reshape_front = 1;
else
    reshape_front = sz_el(3:end);
end

coords = reshape(coords, [reshape_front sz_b(3)]);
