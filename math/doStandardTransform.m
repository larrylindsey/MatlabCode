function y = doStandardTransform(T, x, o, fun, ~)
A = fun(x, [], o);

y = A * T;
end
