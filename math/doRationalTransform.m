function y = doRationalTransform(T, x, o, ~, pfun)

fun = pfun.fun;

if numel(o) == 1
    o = [o o];
end

A_N = fun(x, [], o(1));
A_D = fun(x, [], o(2));
A_D(:,1) = [];

T_N = T(1:size(A_N,2),:);
T_D = T;
T_D(1:size(A_N,2),:) = [];

y_N = A_N * T_N;
y_D = A_D * T_D;

y = y_N ./ (1 + y_D);
end
