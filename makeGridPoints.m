function pts = makeGridPoints(split, tVal)

tol = 2 * pi / 32;

pts = [];

for i_s = 1:numel(split)
   imT = split(i_s).im > tVal;
   [H T R] = hough(imT);
   P = houghpeaks(H, 8);
   lines = houghlines(imT, T, R, P, 'FillGap', 32, 'MinLength', 8);
   pt1 = [lines.point1];
   pt2 = [lines.point2];
   
   pt1 = reshape(pt1, [2 numel(pt1) / 2]);
   pt2 = reshape(pt2, [2 numel(pt2) / 2]);
   
   xy = pt2 - pt1;
   
   angle = atan2(xy(2,:), xy(1,:));
   
   m = pi / 2;
   
   angleModQ = mod(angle + m / 2, m) - m / 2;
   
   lselect = logical(abs(angleModQ) < tol);      
     
   m = pi;
   
   angleModH_horiz = mod(angle + m / 2, m) - m / 2;
   
   angleModH_vert = mod(angle + m / 2 - pi / 2, m) - m / 2;   
   
   hselect = logical(abs(angleModH_horiz) < tol);
   vselect = logical(abs(angleModH_vert) < tol);
   
   hselect = and(hselect, lselect);
   vselect = and(vselect, lselect);
   
   hlines = lines(hselect);
   vlines = lines(vselect);
   
   if ~isempty(hlines) && ~isempty(vlines)
       origin = [split(i_s).r - 1; split(i_s).c - 1];
       for i_h = 1:numel(hlines)
           for i_v = 1:numel(vlines)
               [pt ok] = intersection(hlines(i_h), vlines(i_v));
               if ok
                   pts = cat(2, pts, pt + origin);
               end
           end
       end
   end

end

end