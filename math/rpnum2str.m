function outStr = rpnum2str(num, n)

if isnan(num)
    outStr = '';
    return;
end

if num > 0
    str = num2str(floor(num));
    k = numel(str) + n;
else
    k = n;
end

if n < 1
    outStr = num2str(round(num));
    return;
end

outStr = num2str(num, k);

dotPos = find(outStr == '.');

if isempty(dotPos)
    nz = 0;
    outStr = [outStr '.'];
else
    nz = numel(outStr) - dotPos;
end

outStr = sprintf('%s%s', outStr, repmat('0', [1 n-nz]));

end