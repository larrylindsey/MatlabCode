function c = covariance(x, y)
mu = mean(x(:));
nu = mean(y(:));

xy = x.*y;
c = mean(xy(:)) - mu * nu;