function rcout = affineAlign(rcfrom, rcto)

tr = regressionTransform(rcfrom, rcto, 1, @legendreMat);

rcout = doTransform(rcfrom, tr);

end