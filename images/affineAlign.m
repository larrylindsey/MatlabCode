function rcout = affineAlign(rcfrom, rcto)

tr = fitTransform(rcfrom, rcto, 1, @taylorMat);

rcout = applyTransform(rcfrom, tr);

end
