function E = qmeshKnn(qmesh, k)

E = [];

for i_q = 1:size(qmesh,1)
    dd = dist2(qmesh(i_q,:), qmesh);
    [junk ss] = sort(dd);
    for i_e = 2:(k + 1)
        E = cat(1, E, [i_q ss(i_e)]);
    end
end


end