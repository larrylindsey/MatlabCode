function br = blend_rectangle(ll, ur, control)
% br = blend_rectangle(zc, oc)
%
%

% Creates a "blend rectangle," which is grid that, for example, may be 
% multiplied element-wise with one of a pair of images to be overlapped and
% blended.  It is zero along the two boundaries abutting the corner
% defined by zc, and one along the two boundaries abutting the one defined
% by oc.  It has the property that the differential along those boundaries
% is always zero, except at the corners between oc and zc where such a
% thing is not possible, and it is smooth in any direction internally.
%
% zc - the vector location of the "zero corner"
% oc - the vector location of hte "one corner"
%
% br - the blend rectangle
%
% %Example, to create a 100x100 blend rectangle, and its complement:
% br1 = blend_rectangle([1 1],[100 100]);
% br1(1,1)
% ans =
%      0
% br1(1,50)
% ans =
%      0
% br1(100,100)
% ans =
%      1
% br1(100,50)
% ans =
%      1
% br2 = blend_rectangle([100 100],[1 1]);
% find(br1 + br2 ~= 1)
% ans =
%    Empty matrix: 0-by-1

warning('off', 'MATLAB:divideByZero');

xmin = ll(1);
xmax = ur(1);
ymin = ll(2);
ymax = ur(2);


% Perform a change-of-variables with x and y so that we end up with a
% unit-square domain, with the requested resolution.
x = xmin:xmax;
x = x - xmin;
x = x / max(x);

y = ymin:ymax;
y = y - ymin;
y = y / max(y);

[xx yy] = meshgrid(x, y);

% Make the blend-rectangle.
br = merge_square(xx, yy, control);

warning('on', 'MATLAB:divideByZero');
end
