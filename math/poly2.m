function [PP, Pstruct] = poly2(xx,yy)

if nargin < 1
    xx = linspace(-1,1,1024);
elseif numel(xx) == 1
    xx = linspace(-1,1,xx);
end

if nargin < 2
    if min(size(xx)) == 1
        [xx,yy] = meshgrid(xx,xx);
    else
        error('yy undefined')
    end
else
    if size(xx) ~= size(yy)
        error('xx and yy must have the same dimensions');
    end
end

%Order 0
P00 = ones(size(xx));

%Order 1
P10 = xx;
P01 = yy;

%Order 2
P20 = xx.^2;
P11 = xx .* yy;
P02 = yy.^2;

PP = cat(3, P00, P10, P01, P20, P11, P02);

nn = load('~/code/matlab/data/poly2norm', 'norm');
norm = nn.norm; 

norm = repmat(reshape(norm, [1 1 6]), [size(xx) 1]);

PP = PP ./ norm;

Pstruct = struct('P00', PP(:,:,1),...
    'P10', PP(:,:,2),...
    'P01', PP(:,:,3),...
    'P20', PP(:,:,4),...
    'P11', PP(:,:,5),...
    'P02', PP(:,:,6));