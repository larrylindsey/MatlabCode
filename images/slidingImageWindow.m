function out = slidingImageWindow(im, sz, fun, p)

if nargin < 4
    p = false;
end

imsz = size(im);
grsel = 1:sz(1);
gcsel = 1:sz(2);
grsel = grsel - floor(mean(grsel));
gcsel = gcsel - floor(mean(gcsel));
out = zeros(imsz(1), imsz(2));

rim = imsz(1);
cim = imsz(2);

if p
    for ii = 1:imsz(1)
        rsel = imod(grsel + ii, rim);
        winr = im(rsel, :, :);
        parfor jj = 1:imsz(2)
            csel = imod(gcsel + jj, cim);
            win = winr(:, csel, :);
            out(ii,jj) = fun(win);
        end
    end
else
    for ii = 1:imsz(1)
        rsel = imod(grsel + ii, rim);
        winr = im(rsel, :, :);
        for jj = 1:imsz(2)
            csel = imod(gcsel + jj, cim);
            win = winr(:, csel, :);
            out(ii,jj) = fun(win);
        end
    end
end


end

function n = imod(m, s)

n = mod(m, s);
n(n == 0) = s;

end