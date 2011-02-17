function c = correlation(x, y)
cov = covariance(x,y);
varx = mean(x(:).^2) - mean(x(:)).^2;
vary = mean(y(:).^2) - mean(y(:)).^2;
c = cov / sqrt(varx * vary);