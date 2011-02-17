function dout = distanceTransform1d(bwin)
% dout = distanceTransform1d(bwin)
% Calculates the distance transform of a 1D input BW vector
% 
% bwin - the input vector.  Must be a boolean image.
% dout - the distance transform of bwin, in the same dimensions as bwin.

doutr = doTrans(bwin);
doutl = doTrans(bwin(end:-1:1));
dout = min(doutr, doutl(end:-1:1));

end

function dout = doTrans(bwin)
dout = zeros(size(bwin));

cnt = inf;
for i_b = 1:numel(bwin)
   if bwin(i_b)
       cnt = 0;
   else
       cnt = cnt + 1;
   end
   dout(i_b) = cnt;
end

end