function hues = wheelHue(atans)
% Returns the direction-mapped hue, 
hues = zeros(size(atans));

xi = [-1 -.5 0 .5 1];
hi = [0 1/6 1/3 2/3 1];


hues(:) = interp1(xi, hi, atans(:) / pi, 'linear');


end