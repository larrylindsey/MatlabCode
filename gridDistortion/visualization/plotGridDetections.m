function plotGridDetections(im, rc_found, varargin)

if nargin < 3
    varargin{1} = 'b.';
end

imshow(im, 'XData', [-1 1], 'YData', [-1 1]);
hold on;
plot(rc_found(:,2), rc_found(:,1), varargin{:});
end