function y = applyInverseTransform(x, tr)

if isempty(tr.inv)
    tr = fitInverseTransform(tr);
end

y = applyTransform(x, tr.inv);

end
