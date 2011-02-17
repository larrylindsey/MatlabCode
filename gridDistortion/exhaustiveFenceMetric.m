function [measure vout] = exhaustiveFenceMetric(m_in, L)
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
L_all = 1:max(L(:,2));
L_include = L(v,2);
L_none = logical(L_include == 0);
L_include = L_include(not(L_none));
%nExcluded = sum(L_none);

pk_include = zeros(size(L_include));

for i_ls = 1:numel(L_include)
    sel = logical(L(:,2) == L_include(i_ls));
    pk_include(i_ls) = max(L(sel,3));
end

L_exclude = setdiff(L_all, L_include);
pk_exclude = zeros(size(L_exclude));

for i_lus = 1:numel(L_exclude)
    sel = logical(L(:,2) == L_exclude(i_lus));
    pk_exclude(i_lus) = max(L(sel, 3));
end

vout = var(pk_include);
%measure = (1 + max(pk_include) - min(pk_include)) * (1 + nExcluded) .* (sum(pk_exclude) + 1) ...
%    .* (var(pk_include).^2 + 1) ./ (sum(pk_include) + 1);

%pkRange = max(pk_include) - min(pk_include);

%* (sum(pk_exclude) + 1)
%measure = (1 + pkRange ) * (sum(pk_exclude) + 1) ...
%    .* (var(pk_include).^2 + 1) ./ (mean(pk_include).^2);

measure = (1 + vout + sum(pk_exclude)/sum(pk_include) + sum(L_none)) ./...
    (mean(pk_include).^2);

end