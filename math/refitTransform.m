function tr = refitTransform(fpts, tpts, tr0)
tr = fitTransform(fpts, tpts, tr0.order, tr0.matrixFun, tr0.param);
end
