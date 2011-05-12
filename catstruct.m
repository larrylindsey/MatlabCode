function sout = catstruct(s1, s2)
sout = s1;
sout(end + 1) = s1(end);
ff = fieldnames(sout);

for ii = 1:numel(ff)
   sout(end).(ff{ii}) = s2.(ff{ii});
end
