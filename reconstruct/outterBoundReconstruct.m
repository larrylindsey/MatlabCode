function [x y] = outterBoundReconstruct(instr)

[xx yy] = reconstructDomainBounds(instr);

x(1) = min(xx(1,:));
x(2) = max(xx(2,:));

y(1) = min(yy(1,:));
y(2) = max(yy(2,:));