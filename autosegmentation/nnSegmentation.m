function imseg = nnSegmentation(imin, sub)
global im;

im = imin;

rWiggle = [-1 0 0 1 ];
cWiggle = [0 -1 1 0];

if numel(size(im)) > 2
    im = rgb2gray(im);
end

sz = size(im);
%sz = sz(1:2);
imseg = zeros(sz);

iSeg = 1;

pos = [1 1];

while all(pos <= sz)
    
    Q = pos;
    
    while size(Q, 1) > 0
        currPos = Q(1,:);
        Q = Q(2:end,:);
        
        r = currPos(1);
        c = currPos(2);
        if imseg(r, c) == 0
            imseg(r, c) = iSeg;
            rNbd = r + rWiggle;
            cNbd = c + cWiggle;

            % Reject r's and c's that would cause array out-of-bounds
            sel = rNbd > 0 & rNbd <= sz(1) & cNbd > 0 & cNbd <= sz(2);
            rNbd = rNbd(sel);
            cNbd = cNbd(sel);
            rNbd = rNbd(:);
            cNbd = cNbd(:);
            
            Q = catQ(Q, r, c, rNbd, cNbd, sub);
            
        end
    end
    
    iSeg = iSeg + 1;
    
    pos = incPos(pos, sz);
    while  all(pos <= sz) && imseg(pos(1), pos(2)) > 0
        pos = incPos(pos, sz);
    end
    
%     if pos(2) >= sz(2) - 1
%         keyboard;
%     end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = incPos(pos, sz)
pos(1) = pos(1) + 1;
if pos(1) > sz(1)
    pos(1) = 1;
    pos(2) = pos(2) + 1;
    disp(pos(2))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = catQ(Q, r, c, rNbd, cNbd, sub)
global im;

nbdVal = im(sub2ind(size(im), rNbd, cNbd));
val = repmat(im(r, c), size(nbdVal));

distance = abs(nbdVal - val) - sub;
distance(distance < 0) = 0;
dmin = min(distance);

sel = distance == dmin;

Q = cat(1, Q, cat(2, rNbd(sel(:)), cNbd(sel(:))));

end