function tr = transformToTR(trans)

if isempty(trans.Tinv)
    trans.Tinv = invertTransform(trans.T, type, trans.data.n, ...
        trans.data.u);
end

tr = maketform('custom', trans.iDim, trans.oDim, @doCubicTransform, ...
    @doInverseCubicTransform, trans);