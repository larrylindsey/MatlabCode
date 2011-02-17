function ms = merge_square(x, y, control)
% control in [l u r d]
warning('off', 'MATLAB:divideByZero');
distl = x;
distr = 1-x;
distd = y;
distu = 1-y;

stackset = cat(3, distl, distu, distr, distd);

stackset = stackset(:,:,logical(control >= 0));
control = control(logical(control >= 0));

zeroset = stackset(:,:,not(logical(control)));
oneset = stackset(:,:,logical(control));

zerodist = min(zeroset, [], 3);
onedist = min(oneset, [], 3);

ms = zerodist ./ (onedist + zerodist);

ms(isnan(ms(:))) = .5;

ms = .5 - .5 * cos(pi * ms);


warning('on', 'MATLAB:divideByZero');