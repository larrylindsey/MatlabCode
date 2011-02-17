function [measure vout]= fenceMetric(m_in, samples, L)
% L in [p 3]
% 1st index - tt vector
% 2nd index - label
% 3rd index - weight
%
% samples in [n 2]
% 1st index - tt vector
% 2nd index - label
v = m_in(1):m_in(2):size(L, 1);
v = round(v);

%val_sel = L(v, 1);
if true
    L_sel = 1:max(L(:,2));
    L_samp = L(v,2);
else
    L_uni = L(v, 2);
    L_samp = samples(:,2);
    L_sel = intersect(L_uni, L_samp);
end

pk_sel = zeros(size(L_sel));

for i_ls = 1:numel(L_sel)
    sel = logical(L(:,2) == L_sel(i_ls));
    pk_sel(i_ls) = max(L(sel,3));
end

L_unsel = setdiff(L_sel, L_samp);
pk_unsel = zeros(size(L_unsel));

for i_lus = 1:numel(L_unsel)
    sel = logical(L(:,2) == L_unsel(i_lus));
    pk_unsel(i_lus) = max(L(sel, 3));
end

vout = var(pk_sel);
measure = (1 + sum(L(v,3) == 0)) .* (sum(pk_unsel) + 1) ...
    .* (var(pk_sel) + 1) ./ (sum(pk_sel) + 1);
end

