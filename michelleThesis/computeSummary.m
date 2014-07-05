function stat = computeSummary(stat_m)
global k name_k

stat.afile = {stat_m.afile};
stat.mfile = {stat_m.mfile};
stat.animal = stat_m(1).animal;
stat.totpix = sum([stat_m.totpix]);

for i_k = 1:k
    cname = name_k{i_k};
    c_cat = [stat_m.(cname)];
    stat.(cname).a = [c_cat.a];
    stat.(cname).m = [c_cat.m];
end

end