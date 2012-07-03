function rcout = affineAlign(rcfrom, rcto)

tr = fitTransform(rcfrom, rcto, 1, @legendreMat);

rcout = applyTransform(rcfrom, tr);

end
