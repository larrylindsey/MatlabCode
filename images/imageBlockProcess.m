function [out rstep cstep] = imageBlockProcess(im, bsz, ovlp, fun, template)

if isempty(ovlp)
    ovlp = 0;    
end

if numel(ovlp) == 1
    if ovlp >= 1 || ovlp < 0
        error('Invalid overlap: %g. Valid overlap in [0 1)', ovlp);
    end
    step = round(bsz * (1 - ovlp));
    step(step < 1) = 1;
else
    step = ovlp;
end

imsz = size(im);

rstep = mkstep(step(1), imsz(1));
cstep = mkstep(step(2), imsz(2));

rsel = mksel(bsz(1));
csel = mksel(bsz(2));

if nargin < 5
    out = zeros(numel(rstep), numel(cstep));
else
    out = repmat(template, [numel(rstep), numel(cstep)]);
end

for ii = 1:numel(rstep)
    rr = imod(rsel + rstep(ii), imsz(1));
    mods = imsz(2);
    imrr = im(rr,:);
    parfor jj = 1:numel(cstep)        
        cc = imod(csel + cstep(jj), mods);
        subim = imrr(:, cc); %#ok<PFBNS>
        out(ii,jj, :) = fun(subim); %#ok<PFBNS>
    end
end


end

function s = mkstep(t, sz)
s = 1:t:sz;
% s = s + ceil((sz - s(end)) / 2);
end

function sel = mksel(sz)
sel = 1:sz;
sel = sel - 1;
% sel = sel - floor(mean(sel));
end

function n = imod(m, s)

n = mod(m, s);
n(n == 0) = s;

end
