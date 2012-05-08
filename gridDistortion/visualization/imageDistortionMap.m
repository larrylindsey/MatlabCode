function imageDistortionMap(tr, n, lim, f, varargin)
%  imageDistortionMap(tr, n, lim, f, ...)
% Creates a visual representation of the distortion map of a transform
% Two subplots are created, one to display the pixel-magnitude of the
% distortion map, and another to display the orientation.
%  
% tr [struct] - the transform struct
% n [1] or [2] - n(1) is the number of sample points, in each dimension. 
%       IE, an n-by-n distortion image will be created.
%       n(2) is the number of sample points to use for the vector map, in
%       other words n(2) * n(2) arrows will be plotted over the
%       color-vector-field. When n is a singleton, the value in n(1) is
%       used instead.
% lim [2] - The boundaries of the domain. For instance, the transform is
%       usually calculated over an image that has been geometrically
%       normalized, so that its pixels are located between -1 and 1, in
%       both x- and y-directions. In this case, lim = [-1 1].
% f [1] - the size factor. In the example above, if the original image was
%       a 16384-by-16384 image, then f = 16384 / (1 - (-1)) = 16384 / 2 =
%       8192.
% ... - any additional arguments are passed directly into the quiver
%       command, used to draw the arrows on the plot


if isempty(lim)
    lim = [-1 1];
end

if numel(n) == 1
    n = [n n];
end

[thresh rcPoly rcFound varargin] = getThreshArg(varargin);

[rcMap, rctrMap] = calculateRC(tr, n(1), lim);

if isempty(rcFound)
    [rcQ, rctrQ] = calculateRC(tr, n(2), lim);
    dfmag0 = [];
    dfcolor0 = [];
    mapSel = [];
else
    mapSel = inpolygon(rcMap(:,1), rcMap(:,2), rcPoly(:,1), rcPoly(:,2));
    
    [dfmag0 dfcolor0] = makeDistortionMap(rcMap, rctrMap, ...
        [], [n(1) n(1)]);
    
    rcMap(~mapSel,:) = nan;
    rctrMap(~mapSel,:) = nan;
    
    rcQ = rcFound;
    rctrQ = doTransform(rcQ, tr);
end

if ~isempty(thresh)
    threshf =  thresh / f;
    tlim = [0 threshf];
else 
    threshf = [];
    tlim = [];
end

[dfmag dfcolor] = makeDistortionMap(rcMap, rctrMap, threshf, [n(1) n(1)]);


figure;
imagesc(dfmag * f, 'XData', lim, 'YData', lim); axis image; axis xy;
colormap(jet(256));
if ~isempty(thresh)
    caxis([0 thresh]);
end
colorbar;

if ~isempty(mapSel)
    mapSel3 = repmat(mapSel, [1 1 3]);
    dfcolor(not(mapSel3)) = dfcolor0(not(mapSel3)) / 2;
    
    dfmag0 = applyColorMap(dfmag0, jet(256));
    dfmag = applyColorMap(dfmag, jet(256), tlim);
    
    dfmag(not(mapSel3)) = dfmag0(not(mapSel3)) / 2;
    hold on;
    image(dfmag, 'XData', lim, 'YData', lim);
end



figure;
imagesc(dfcolor, 'XData', lim, 'YData', lim); axis image; axis xy;
hold on;
quiver(rcQ(:,2), rcQ(:,1), ...
    rcQ(:,2) - rctrQ(:,2), rcQ(:,1) - rctrQ(:,1), varargin{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RC RCtr] = calculateRC(trans, n, lim)
lilr = linspace(lim(1), lim(2), n);
lilc = lilr;

[R C] = meshgrid(lilr, lilc);
RC = cat(2, C(:), R(:));
RCtr = doTransform(RC, trans);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t rcHull rc incell] = getThreshArg(incell)

rcHull = [];

[argstr incell] = getargs({'rc', 'threshold'}, incell);

rc = argstr.rc;
t = argstr.threshold;

if  ~isempty(rc)
    k = convhull(rc(:,2), rc(:,1));
    rcHull = rc(k,:);
end

end
