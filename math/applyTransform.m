function y = applyTransform(x, tr)

y = tr.doTrans(tr.T, x, tr.order, tr.matrixFun, tr.param);

end
