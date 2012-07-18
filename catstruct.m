function sout = catstruct(s1, s2)
if numel(s1) < 1
    sout = s2;
else
    
    sout = s1;
    sout(end + 1) = s1(end);
    ff = fieldnames(sout);
    
    for ii = 1:numel(ff)
        sout(end).(ff{ii}) = s2.(ff{ii});
    end
end
