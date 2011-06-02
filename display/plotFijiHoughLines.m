function [x, y] = plotFijiHoughLines(r, t, im, nr, nt, ts)
figure;
imshow(im);
hold on;

diag = sqrt(sum(size(im).^2));
ll = floor(-diag):ceil(diag);
nl = numel(ll);
ll = repmat(ll, [2, 1]);

if nargin < 6
    ts = 0;
    if nargin < 5
        nt = [];
        if nargin < 4
            nr = [];
        end
    end       
end

if isempty(nt)
    nt = 180;
end

if isempty(nr)
    nr = 2 * floor(diag);
end

rho = linspace(-diag, diag, nr);
theta = linspace(-pi / 2, pi / 2, nt + 1);
t = t + 1;
r = r + 1;

for ii = 1:numel(r)
    t0 = rho(r(ii)) * [cos(theta(t(ii))); sin(theta(t(ii)))];
    t0 = repmat(t0, [1, nl]);
    tparm = [-sin(theta(t(ii))); cos(theta(t(ii)))];
    tparm = repmat(tparm, [1 nl]);
    
    disp(180 * theta(t(ii)) / pi);
    disp(rho(r(ii)));
    
    linepath = tparm .* ll + t0;
    
    x = linepath(1,:);
    y = linepath(2,:);
    
    plot(x, y);    
        
    if t > 0
       pause(ts); 
    end
        
end


end