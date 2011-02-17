function [x y] = innerBoundReconstruct(instr)

[xx yy] = reconstructDomainBounds(instr);

x(1) = max(xx(1,:));
x(2) = min(xx(2,:));

y(1) = max(yy(1,:));
y(2) = min(yy(2,:));