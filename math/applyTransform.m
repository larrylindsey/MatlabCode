function y = applyTransform(x, tr)
if isempty(x)
    y = [];
else
    y = tr.doTrans(tr.T, x, tr.order, tr.matrixFun, tr.param);
end

end
