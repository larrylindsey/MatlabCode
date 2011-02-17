function dout = sqFloatDistanceTransform1d(bwin, fvect)

doutr = doTrans(bwin, fvect);
doutl = doTrans(bwin(end:-1:1), fvect(end:-1:1));
dout = min(doutr, doutl(end:-1:1));

end

function dout = doTrans(bwin, fvect)
dout = zeros(size(bwin));

cnt = inf;
s0 = inf;
for i_b = 1:numel(bwin)
   if bwin(i_b) ~= 0
       cnt = 0;
       s0 = fvect(i_b);
   else
       cnt = cnt + 1;
   end
   dout(i_b) = s0 + cnt;
end

end