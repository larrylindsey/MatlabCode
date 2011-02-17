function val = trapint2(xx,yy,ff)
% val = trapint2(xx, yy, ff)
%
% Calculates a (rather error prone) numerical trapezoidal integral of
% the function ff(xx, yy)

%Make sure that all matrices have the same size
if size(xx) ~= size(yy)
    error('size(xx) ~= size(yy)');
elseif size(yy) ~= size(ff)
    error('size(yy) ~= size(ff)');
end

%Trapezoidal integration weights.
ww = ones(size(xx));
ww([1 end],:) = .5;
ww(:,[1 end]) = .5;
ww([1 end], [1 end]) = .25;

%Get the trapezoid widths
dx = diff(xx, 1, 2);
dy = diff(yy, 1, 1);

%We will interpolate dx and dy so that they are a function of xx, yy.

%Here, find the midpoints of the x-values in the x-direction, and
%y-values in the y-direction.  These values, mx, my, will be the original
%domains for interpolating dx, dy.  Ie, we want dx(mx) -> dx(xx), and
%dy(yy) -> dy(my)
mx = .5 * (xx(:, 1:(end-1)) + xx(:,2:end));
my = .5 * (yy(1:(end-1),:) + yy(2:end,:));

%Interpolate!
dx = interp2(mx, yy(:,2:end), dx, xx, yy,'linear');
dy = interp2(xx(2:end,:), my, dy, xx, yy,'linear');
%interp2 doesn't do extrapolation, instead replacing extrapolated values
%with Nan.  Here, we integrate over actual numerical values (ie, non-nans),
%then scale up to do a real crappy extrapolation.

val = ff .* ww .* dx .* dy;
num_select = not(isnan(val)); %Select non-nans
num_count = sum(num_select(:)); %Count non-nans
val = val * numel(val) / num_count; %Scale by #values / #non-nans

val = sum(val(num_select(:)));
