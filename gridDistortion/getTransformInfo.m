function [fwdFunc, invFunc, T, transinfo] = getTransformInfo(r, c, rTest, cTest)

r = condition(r);
c = condition(c);
rTest = condition(rTest);
cTest = condition(cTest);

universe = cat(2, r, c, rTest, cTest);

[outSpace junk emin] = ransac(universe, @makeCubicModel, ...
    @measureModelCubic, .1, ceil(size(universe, 1) / 4), 12, 1000); %#ok

%keyboard

%pr = leastSquaresQuadraticFit2(r, c, dv(:, 1));
%pc = leastSquaresQuadraticFit2(r, c, dv(:, 2));

%Tf = [pr pc];
%Tf(4:5, :) = Tf(4:5, :) + eye(2);

ransacR = outSpace(:, 1)';
ransacC = outSpace(:, 2)';

ransacRTest = outSpace(:, 3);
ransacCTest = outSpace(:, 4);

pr = leastSquaresONCubicFit2(ransacR, ransacC, ransacRTest);
pc = leastSquaresONCubicFit2(ransacR, ransacC, ransacCTest);

Ti = [pr pc];

transinfo.emin = emin;
transinfo.r = ransacR;
transinfo.c = ransacC;

pr = leastSquaresONCubicFit2(ransacRTest, ransacCTest, ransacR);
pc = leastSquaresONCubicFit2(ransacRTest, ransacCTest, ransacC);

Tf = [pr pc];
%Ti(4:5, :) = Ti(4:5, :) + eye(2);

T = cat(1, Tf, Ti);

invFunc = doCubicTransform2(false, 'legendre');
fwdFunc = doCubicTransform2(true, 'legendre');

end
% 
% function pp = makeQuadraticModel(in, junk)
% pr = leastSquaresQuadraticFit2(in(:,1)', in(:,2)', in(:,3));
% pc = leastSquaresQuadraticFit2(in(:,1)', in(:,2)', in(:,4));
% pp = [pr pc];
% %pp(4:5, :) = pp(4:5, :) + eye(2);
% end

function pp = makeCubicModel(in, junk)
pr = leastSquaresONCubicFit2(in(:,1)', in(:,2)', in(:,3));
pc = leastSquaresONCubicFit2(in(:,1)', in(:,2)', in(:,4));
pp = [pr pc];
%pp(4:5, :) = pp(4:5, :) + eye(2);
end
% 
% function e = measureModel(pp, in, junk)
% U = doQuadraticTransform2(in(:,1:2), pp);
% dd = (U - in(:,3:4)) ./ repmat(sqrt(sum(in(:,1:2).^2, 2)), [1 2]);
% e = sqrt(mean(dd(:).^2));
% end

function e = measureModelCubic(pp, in, junk)
U = doCubicTransform2(in(:,1:2), pp, 'legendre');
dd = (U - in(:,3:4)) ./ repmat(sqrt(sum(in(:,1:2).^2, 2)), [1 2]);
e = sqrt(mean(dd(:).^2));
end

function in = condition(in)

if size(in, 1) < size(in, 2)
    in = in';
end

end