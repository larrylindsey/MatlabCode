function plotRCComparison(rc1, rc2, argsAll, args1, args2)

if nargin < 5
    args2 = 'kx';
    if nargin < 4
        args1 = 'bx';
        if nargin < 3
            argsAll = {};
        end
    end
end

if ~iscell(args2) 
    args2 = {args2};
end
if ~iscell(args1)
    args1 = {args1};
end
if ~iscell(argsAll)
    argsAll = {argsAll};
end

rc1 = shapeRC(rc1);
rc2 = shapeRC(rc2);

r = cat(2, rc1(:,1), rc2(:,1));
c = cat(2, rc1(:,2), rc2(:,2));

plot(r', c', argsAll{:});
hold on;
if ~isempty(args1{1})
    plot(rc1(:,1), rc1(:,2), args1{:});
end
if ~isempty(args2{1})
    plot(rc2(:,1), rc2(:,2), args2{:});
end

grid on;
axis image;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rc = shapeRC(rc)
if size(rc, 2) > size(rc, 1)
    rc = rc.';
end
end