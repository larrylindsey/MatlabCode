function [rlReg clReg] = findRegularCrossPoints(r, c, rl, cl)

rl = normlabel(rl);
cl = normlabel(cl);

rows = cell(1, max(rl(:)));
cols = cell(1, max(cl(:)));
rlReg = zeros(1, max(rl(:)));
clReg = zeros(1, max(cl(:)));

for i_x = 1:numel(c)
    i_r = rl(round(r(i_x)));
    i_c = cl(round(c(i_x)));
    
    if i_r == 0
        warning('findRegularCrossPoints:mislabel', ...
            ['Found a cross point that does not appear to fall within a'...
            ' legitimate label']);
    else
        rows{i_r} = r(i_x);
    end
    
    if i_c == 0
        if i_r ~= 0 %No reason to send the same warning twice.
            warning('findRegularCrossPoints:mislabel', ...
                ['Found a cross point that does not appear to fall '...
                'within a legitimate label']);
        end
    else
        cols{i_c} = c(i_x);
    end        
end

for i_r = 1:numel(rlReg)
    rlReg(i_r) = median(rows{i_r});
end

for i_c = 1:numel(clReg)
    clReg(i_c) = median(cols{i_c});
end

rlReg = linefit(rlReg);
clReg = linefit(clReg);


end

function lout = normlabel(lin)

cnt = 1;

lout = lin;

for i_l = 1:max(lin)
    sel = logical(lin == i_l);
    if any(sel)
        lout(sel) = cnt;
        cnt = cnt + 1;
    end
end

end


function fit = linefit(unfit)

x = 1:numel(unfit);
sz = size(unfit);
if sz(2) > sz(1)
    unfit = unfit';
end

A = cat(2, x', ones(size(x')));

AtA = A' * A;
[U S V] = svd(AtA);

Ssel = logical(eye(size(S)));
Sinv = zeros(size(S));
Sinv(Ssel) = 1./ (S(Ssel));
AtAinv = U * Sinv * V';

c = AtAinv * (A') * unfit;

fit = c(1) * x + c(2);

end