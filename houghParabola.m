function [h adom bdom cdom]= houghParabola(imbw, ap, bp, cp)

[yr xc] = find(imbw);

alim = ap(1);
ares = ap(2);
%ga = 1 / ares;

blim = bp(1);
bres = bp(2);
%gb = 1 / bres;

clim = cp(1);
cres = cp(2);
%gc = 1 / cres;



adom = linspace(-alim, alim, ares);
bdom = linspace(-blim, blim, bres);
cdom = linspace(0, clim, cres);


h = zeros(ares, cres);

%disp(numel(xc));

sprintf('init');
lastN = 4;

for i_x = 1:round(numel(xc) / 4):numel(xc)
    x = xc(i_x);
    y = yr(i_x);
    a = (y - cdom) / (x * x);
    
    displaystr = sprintf('%g x of %g - %2.1f %%', i_x, numel(xc), i_x / numel(xc));
    rmstr = repmat('\b', [1, lastN]);
    lastN = numel(displaystr);
    fprintf(rmstr);
    fprintf(displaystr);
    
    
    i_c = find(and(logical(a >= -alim), logical(a <= alim)));
    
    for ii = i_c
        i_a = find(a(ii) < adom,1, 'first');
        h(i_a, i_c) = h(i_a, i_c) + 1;
    end
end

end
