function trOut = constrainTransform(trIn)
% Right now, works for taylorMat transforms

rotProj =...
   [0   0;      % Order 0
    0   3/8;    % Order 1
    3/8 0;
    0   0;      % Order 2
    0   0;
    0   0;
    0   3/4;    % Order 3
    3/4 0;
    0   3/4;
    3/4 0;
    0   0;      % Order 4
    0   0;
    0   0;
    0   0;
    0   0;
    0   -1/8;   % Order 5
    -1/8    0;
    0   -1/4;
    -1/4    0;
    0   -1/8;
    -1/8    0];

idProj = zeros(size(rotProj));
idProj(1:3,:) = ...
   [0   0;
    1   0;
    0   1];
    

pbProj = zeros(size(rotProj));
pbProj(7, 1) = 1;
pbProj(10,2) = 1;

T = trIn.T;

trOut = trIn;

Trot = doProj(T, fliplr(rotProj));
Tid = doProj(T, fliplr(idProj));
Tpb = doProj(T, fliplr(pbProj));

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
tr = fitTransform(xyaff, xyt, tr.order, @legendreMat, tr.data);
end
