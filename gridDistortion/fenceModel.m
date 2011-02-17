function mout = fenceModel(samples, junk)
sx = samples(:,1);

sx = sort(sx);
d = min(diff(sx));
s = min(sx);

s = mod(s, d);
if s == 0
    s = d;
end

mout = [s d];

end