function trOut = constrainTransform2(trIn)
% Right now, works for taylorMat transforms
origMorder = [...
    0 0;
    1 0;
    0 1;
    2 0;
    1 1;
    0 2;
    3 0;
    2 1;
    1 2;
    0 3;
    4 0; 
    3 1;
    2 2;
    1 3;
    0 4;
    5 0;
    4 1;
    3 2;
    2 3;
    1 4;
    0 5];

rotProj =...
   [0   0;      % Order 0
    0   0;    % Order 1
    0   0;
    0   0;      % Order 2
    0   0;
    0   0;
    0   1;    % Order 3
    1   0;
    0   1;
    1   0;
    0   0;      % Order 4
    0   0;
    0   0;
    0   0;
    0   0;
    0   0;   % Order 5
    0   0;
    0   0;
    0   0;
    0   0;
    0   0];

idProj = zeros(size(rotProj));
idProj(1:3,:) = ...
   [0   0;
    1   0;
    0   1];
    

pbProj = zeros(size(rotProj));
pbProj(7, 1) = 1;
pbProj(10,2) = 1;

pbProj = rectifyOrder(pbProj, origMorder, 5);
idProj = rectifyOrder(idProj, origMorder, 5);
rotProj = rectifyOrder(rotProj, origMorder, 5);

T = trIn.T;

trOut = trIn;

Trot = doProj(T, rotProj);
Tid = doProj(T, idProj);
Tpb = doProj(T, pbProj);

trOut.T = Trot + Tid + Tpb;
%trOut.Tinv = [];

trOut = affinize(trOut);

% trOut = populateTransInverse(trOut);

end

function Tout = doProj(T, p)
p = p(1:size(T,1),:);

k = dot(T(:),p(:)) / dot(p(:),p(:));

Tout = k * p;
end

function tr = affinize(tr)
n = 32;
[X Y] = meshgrid(linspace(tr.data.u(1), tr.data.u(2), n),...
    linspace(tr.data.v(1), tr.data.v(2), n));
xy = cat(2, X(:), Y(:));
xyt = applyTransform(xy, tr);
xyaff = affineAlign(xy, xyt);
tr = fitTransform(xyaff, xyt, tr.order, @legendreMat);
end
