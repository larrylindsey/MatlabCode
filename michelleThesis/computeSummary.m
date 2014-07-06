function stat = computeSummary(stat_m)


stat.afile = {stat_m.afile};
stat.mfile = {stat_m.mfile};
stat.animal = stat_m(1).animal;
stat.totpix = sum([stat_m.totpix]);

fields = fieldnames(stat_m);

for i_f = 1:numel(fields)
    cname = fields{i_f};
    if isstruct([stat_m.(cname)])
        c_cat = [stat_m.(cname)];
        stat.(cname).a = [c_cat.a];
        stat.(cname).m = [c_cat.m];
    end
end

end