function [tr1 tr2] = constrainTransformIxy(trIn)
% Right now, works for taylorMat transforms

tr1 = projectD1(trIn);
tr2 = projectD2(trIn);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = projectD1(trIn)
rotProjX =...
   [0   0;      % Order 0
    0   0;    % Order 1
    3/8 0;
    0   0;      % Order 2
    0   0;
    0   0;
    0   0;    % Order 3
    3/4 0;
    0   0;
    3/4 0;
    0   0;      % Order 4
    0   0;
    0   0;
    0   0;
    0   0;
    0   0;   % Order 5
    -1/8    0;
    0   0;
    -1/4    0;
    0   0;
    -1/8    0];

rotProjY =...
   [0   0;      % Order 0
    0   3/8;    % Order 1
    0   0;
    0   0;      % Order 2
    0   0;
    0   0;
    0   3/4;    % Order 3
    0   0;
    0   3/4;
    0   0;
    0   0;      % Order 4
    0   0;
    0   0;
    0   0;
    0   0;
    0   -1/8;   % Order 5
    0   0;
    0   -1/4;
    0   0;
    0   -1/8;
    0   0];


idProj = zeros(size(rotProjX));
idProj(1:3,:) = ...
   [0   0;
    1   0;
    0   1];
    

pbProjX = zeros(size(rotProjX));
pbProjY = pbProjX;
pbProjX(7, 1) = 1;
pbProjY(10,2) = 1;

T = trIn.T;

tr = trIn;

TrotX = doProj(T, rotProjX);
TrotY = doProj(T, rotProjY);
Tid = doProj(T, idProj);
TpbX = doProj(T, pbProjX);
TpbY = doProj(T, pbProjY);

tr.T = TrotX + TrotY + Tid + TpbX + TpbY;


tr = affinize(tr);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = projectD2(trIn)
rotProjX = zeros(21,2);
rotProjY = rotProjX;

rotProjX([8, 10],1) = 1;
rotProjY([7, 9], 2) = 1;

idProj = zeros(size(rotProjX));
idProj(1:3,:) = ...
   [0   0;
    1   0;
    0   1];
    

pbProjX = zeros(size(rotProjX));
pbProjY = pbProjX;
pbProjX(7, 1) = 1;
pbProjY(10,2) = 1;

T = trIn.T;

tr = trIn;

TrotX = doProj(T, rotProjX);
TrotY = doProj(T, rotProjY);
Tid = doProj(T, idProj);
TpbX = doProj(T, pbProjX);
TpbY = doProj(T, pbProjY);

tr.T = TrotX + TrotY + Tid + TpbX + TpbY;


tr = affinize(tr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tout = doProj(T, p)
p = p(1:size(T,1),:);

k = dot(T(:),p(:)) / dot(p(:),p(:));

Tout = k * p;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = affinize(tr)
n = 32;
[X Y] = meshgrid(linspace(tr.data.u(1), tr.data.u(2), n),...
    linspace(tr.data.v(1), tr.data.v(2), n));
xy = cat(2, X(:), Y(:));
xyt = doTransform(xy, tr);
xyaff = affineAlign(xy, xyt);
tr = regressionTransform(xyaff, xyt, tr.order, @legendreMat, tr.data);
end
