function strout = findByField(strin, val, f)

if ischar(val)
    eqTest = @strTest;
elseif isnumeric(val)
    eqTest = @numTest;
else
    error('Could not figure out what to do.  Poop.');
end

for i_str = 1:numel(strin)
    str = strin(i_str);
    if eqTest(val, str.(f))
        strout = str;
        return;
    end
end

strout = [];

end

function test = strTest(one, two)
test = ischar(two) & strcmp(one, two);
end

function test = numTest(one, two)
test = isnumeric(two) & one == two;
end