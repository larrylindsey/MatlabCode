function [PP, Pstruct] = orthonormpoly2(xx,yy)

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
P00 = ones(size(xx)) / 2;

%Order 1
P10 = sqrt(.75) * xx;
P01 = sqrt(.75) * yy;

%Order 2
P20 = sqrt(1.25) * .5*(3 * xx.^2 - 1);
P11 = 1.5 * xx .* yy;
P02 = sqrt(1.25) * .5*(3 * yy.^2 - 1);

PP = cat(3, P00, P10, P01, P20, P11, P02);
Pstruct = struct('P00', P00,...
    'P10', P10,...
    'P01', P01,...
    'P20', P20,...
    'P11', P11,...
    'P02', P02);