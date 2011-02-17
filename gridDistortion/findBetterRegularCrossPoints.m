function [rlReg clReg] = findBetterRegularCrossPoints(r, c, rmodel, cmodel)

rlReg = makeItHappen(max(r), rmodel);
clReg = makeItHappen(max(c), cmodel);

end

function out = makeItHappen(xmax, xmodel)
xmodulus = xmodel(2);
xstart = mod(xmodel(1), xmodulus);
xfinish = xmax + xmodulus - mod(xmax, xmodulus);
out = xstart:xmodulus:xfinish;
end
