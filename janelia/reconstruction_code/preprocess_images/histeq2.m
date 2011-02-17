function [out,T] = histeq2(a,image)
% HISTEQ Enhance contrast using histogram equalization.
%
% image = histeq2(image(image>0.1 & image<0.9), image)
%

iptchecknargin(1,3,nargin,mfilename);

NPTS = 256;

n = 64; % Default n
hgram = ones(1,n)*(numel(a)/n);
n = NPTS;
kind = 1;

classChanged = false;
if isa(a,'int16')
  classChanged = true;
  a = im2uint16(a);
end

if min(size(hgram))>1
   msg = 'HGRAM must be a vector.';
   eid = sprintf('Images:%s:hgramMustBeAVector',mfilename);
   error(eid, msg);
end

% Normalize hgram
hgram = hgram*(numel(a)/sum(hgram));       % Set sum = numel(a)
m = length(hgram);

% Compute cumulative histograms
if kind==1,
   nn = imhist(a,n)';
   cum = cumsum(nn);
else
   % Convert image to equivalent gray image
   I = ind2gray(a,cm);
   nn = imhist(I,n)';
   cum = cumsum(nn);
end
cumd = cumsum(hgram*numel(a)/sum(hgram));

% Create transformation to an intensity image by minimizing the error
% between desired and actual cumulative histogram.
tol = ones(m,1)*min([nn(1:n-1),0;0,nn(2:n)])/2;
err = (cumd(:)*ones(1,n)-ones(m,1)*cum(:)')+tol;
d = find(err < -numel(a)*sqrt(eps));
if ~isempty(d)
   err(d) = numel(a)*ones(size(d));
end
[dum,T] = min(err); %#ok
T = (T-1)/(m-1);

if kind == 1 % Modify intensity image
  image_1 = round(image*255);
   b = T(image_1+1);
else % Modify colormap by extending the (r,g,b) vectors.

   % Compute equivalent colormap luminance
   ntsc = rgb2ntsc(cm);

   % Map to new luminance using T, store in 2nd column of ntsc.
   ntsc(:,2) = T(floor(ntsc(:,1)*(n-1))+1)';

   % Scale (r,g,b) vectors by relative luminance change
   map = cm.*((ntsc(:,2)./max(ntsc(:,1),eps))*ones(1,3));

   % Clip the (r,g,b) vectors to the unit color cube
   map = map ./ (max(max(map')',1)*ones(1,3));
end

if nargout==0,
   if kind==1
      imshow(b);
   else
      imshow(a,map),
   end
   return
end

if kind==1
  if classChanged
   out = im2int16(b);
  else
    out = b;
  end
else
   out = map;
end
