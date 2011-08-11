function ph = display_seam(seam, ah, varargin)
%ph = display_seam(seam, ah, <plot args>)
%
%Overlays a seam onto the current image.
%
% seam - the seam to overlay.  For a horizontal seam, this should be a row
%        vector, otherwise a column vector.
% ah - the handle to the axes over which to plot the seam.  If left
%      unspecified, or empty, the value returned by gca will be used.
% Further arguments are passed directly into plot.
%
% ph - the handle to the plot object, as returned by plot
%
% See also plot

if nargin < 2 || isempty(ah)
    ah = gca;
end

%Don't want to erase the image!
hold(ah,'on');

if size(seam,1) > size(seam,2)
    %Seam is a vertical seam
    ph = plot(ah, seam(:), 1:numel(seam), varargin{:});
else
    %Seam is a horizontal seam
    ph = plot(ah, 1:numel(seam), seam(:), varargin{:});
end
    