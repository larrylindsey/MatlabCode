function [bestR,bestP,bestF] = maxF_2(R,P)
% cut out of boundaryBench.m
% interpolate to find best F and coordinates thereof
bestR = R(1);
bestP = P(1);
bestF = fmeasure(R(1),P(1));
for i = 2:numel(R),
  for d = linspace(0,1),
    r = R(i)*d + R(i-1)*(1-d);
    p = P(i)*d + P(i-1)*(1-d);
    f = fmeasure(r,p);
    if f > bestF,
      bestR = r;
      bestP = p;
      bestF = f;
    end
  end
end

% compute f-measure fromm recall and precision
function [f] = fmeasure(r,p)
f = 2*p.*r./(p+r+((p+r)==0));
