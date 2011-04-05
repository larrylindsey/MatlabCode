function [segout ws_imseg ws_imim] = shatterSegment(segin, varargin)
% Performs a seeded watershed segmentation, based on ilastik segment
% prediction and the raw image.
%
% function [segout ws_imseg ws_imim] = better_segment(segin, im)
% segin - the segment prediction image for cytoplasm (membrane should be
%           darker)
% im - the raw em image.
%
% or
% function [segout ws_imseg ws_imim] = better_segment(ilp, [ip])
% ilp - a .ilp or .h5 file corresponding to an ilastik project
% ip - label index corresponding to cytoplasm.
%
%
%

% Author: Larry Lindsey, July 2010

%Default Values
dns = 4;
dimth = 200;
dimbwth = .35;
dsegth = 200;

pth = .4;

if ischar(segin)
    info = hdf5info(segin);
    
    if nargin > 1
        dsel = varargin{2};
        varargin = {varargin{3:end}};
    else
        dsel = 1;
    end
    
    im = hdf5read(info.GroupHierarchy.Groups.Datasets(1));
    im = uint8(squeeze(im));
    
    segin = hdf5read(info.GroupHierarchy.Groups.Datasets(3));
    segin = squeeze(segin);
    segin = squeeze(segin(dsel, :, :));
else
    im = varargin{1};
    varargin = {varargin{2:end}};
end

if ~isempty(varargin)
    ns = varargin{1};
    varargin = {varargin{2:end}};
else  
    ns = dns;
end

if ~isempty(varargin)
    segth = varargin{1};
    varargin = {varargin{2:end}};
else  
    segth = dsegth;
end

if ~isempty(varargin)
    imbwth = varargin{1};
    varargin = {varargin{2:end}};
else  
    imbwth = dimbwth;
end

if ~isempty(varargin)
    imth = varargin{1};
else  
    imth = dimth;
end



if size(segin, 3) > 1
    segin = segin(:,:,1);
end

if ns == 8
    shifts = [1  1;
             -1  1;
             -1 -1;
              1 -1;
              0 1;
              1 0;
              0 -1;
              -1 0];
elseif ns == 4
    shifts = [1  1;
             -1  1;
             -1 -1;
              1 -1];
else
    error('%d invalid value for ns', ns);
end

invSegin = 1 - im2double(segin);
segbw = invSegin < pth;
seginMask = imerode(segbw, strel('disk', 2));
seginMask = bwareaopen(seginMask, segth, 4);


invIm = 1 - im2double(im);
invImMed = medfilt2(invIm, [3 3]);
imMask = bwareaopen(invImMed < imbwth, imth, 4);
im_minim = imimposemin(invImMed, imMask, 4);
ws_imim = watershed(im_minim , 4);

im_minseg = imimposemin(invIm, seginMask, 4);
ws_imseg = watershed(im_minseg, 4);

for is = 1:size(shifts, 1)
     ws_imseg = ws_imseg .* shiftseed(im_minseg, seginMask, ...
         shifts(is, :));
end

%segout = segment_and(ws_imseg, ws_imim, seginMask, 1000);
segout = ws_imseg .* ws_imim;
segout = segout ~= 0;
segout = bwlabel(segout, 4);

end

function wsout = shiftseed(im, mask, sdir)

rsel = 1:size(mask,1);
csel = 1:size(mask,2);
rsel = modi(rsel + sdir(1), rsel(end));
csel = modi(csel + sdir(2), csel(end));

mask = mask(rsel, csel);
im_min = imimposemin(im, mask, 4);

wsout = watershed(im_min, 4);

end