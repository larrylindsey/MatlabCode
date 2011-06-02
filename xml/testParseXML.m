function testParseXML(sel)

test{1} = {'<foo bar="baz"/>'};
test{2} = {'<foo bar="baz">baz="holla"</foo>'};
test{3} = {'<foo', 'bar="baz"', '/>'};
test{4} = {'<foo', 'bar="baz"', '>', 'baz="holla"', '</foo>'};
test{5} = {'<foo> baz="holla"</foo>'};
test{6} = {'<foo>', 'baz="holla"', '</foo>'};

n = numel(test);

if nargin < 1
    sel = 1:n;
end

for i_n = sel
    count = 0;
    [xmst pos name] = parseXML(test{i_n});
    
    while sum(pos > 0)
        count = count + 1;
        [xmst pos name] = parseXML(test{i_n}, pos);
    end
    disp(xmst);
    disp(name);
end


end