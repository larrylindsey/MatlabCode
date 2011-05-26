function [gridRC, model] = fitGrid(model, rc, modelErr)

modelErr = modelErr * max(rms(model, 2)) * sqrt(2);

% Right reference, left reference
lref = model(1,:);
rref = model(2,:);
% refs in [1 2]

% index into rc for the origin-corner-most grid point, with respect to
% lref, rref
[i_corner l_flip r_flip] = findCorner(rc, lref, rref);

if l_flip
    lref = -lref;
end
if r_flip
    rref = -rref;
end


% left shift, left link
[lShift, lLink] = shiftLink(rc, lref, modelErr);

% right shift, right link
[rShift, rLink] = shiftLink(rc, rref, modelErr);


% condition shifts.  For gridRC, the rows correspond to the lref direction,
% columns to rref direction.

% TODO:  Check that orientation is preserved.
rc_lShift = cat(2, lShift, zeros(size(lShift)));
rc_rShift = cat(2, zeros(size(rShift)), rShift);

lsel = lLink ~= 0;
lLink = cat(2, lLink, (1:numel(lLink))');
lLink = lLink(lsel,:);
rc_lShift = rc_lShift(lsel,:);

rsel = rLink ~= 0;
rLink = cat(2, rLink, (1:numel(rLink))');
rLink = rLink(rsel, :);
rc_rShift = rc_rShift(rsel, :);


gridRC = linkGrid(rc, cat(1, rc_lShift, rc_rShift), ...
    cat(1, lLink, rLink), ...
    i_corner);

model(1,:) = lref;
model(2,:) = rref;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gridRC] = linkGrid(rc, shift, link, i_corner)


adjmat = sparse([link(:,1) ; link(:,2)], [link(:,2) ; link(:,1)], ...
    ones(size(link, 1) * 2, 1));

lshift = sparse([link(:,1) ; link(:,2)], [link(:,2) ; link(:,1)], ...
    [-shift(:,1) ; shift(:,1)]);
rshift = sparse([link(:,1) ; link(:,2)], [link(:,2) ; link(:,1)], ...
    [-shift(:,2) ; shift(:,2)]);

gridRC = nan(size(rc));

gridRC(i_corner, :) = 0;

fifo = i_corner;

while ~isempty(fifo)
    
    index = fifo(1);
    fifo = fifo(2:end);
    
    neighbors = find(adjmat(:, index) > 0);
    
    adjmat(:,index) = 0; %#ok<SPRIX>
    adjmat(index,:) = 0; %#ok<SPRIX>
    
    if ~isempty(neighbors)
        offset = cat(2, lshift(neighbors, index), rshift(neighbors, index));
        ctr = repmat(gridRC(index,:), [numel(neighbors) 1]);

        gridRC(neighbors, :) = offset + ctr;

        neighbors = setdiff(neighbors, fifo);

        if isempty(fifo)
            fifo = neighbors;
        else
            fifo = cat(1, fifo, neighbors);
        end
    end
end

gridRC = gridRC - repmat(min(gridRC, [], 1), [size(gridRC,1), 1]);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shift, link] = shiftLink(rc, ref, err)

n = 1;

%rcRefFactor = rms(rc ./ repmat(ref, [size(rc, 1), 1]), 2) * sqrt(2);
%maxShift = ceil(max(rcRefFactor) - min(rcRefFactor));
maxShift = round(rms(ref) * sqrt(2) / err / 3);

shift = zeros(size(rc, 1), 1);
link = shift;

while any(link == 0) && n <= maxShift
    rcShift = rc + repmat(ref * n, [size(rc, 1), 1]);
    dmat = dist2(rc, rcShift);
    [d di] = min(dmat, [], 2);
    sel = and(d < n * (err * err), link == 0);
    link(sel) = di(sel);
    shift(sel) = n;
    n = n + 1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i_corner l_flip r_flip] = findCorner(rc, lref, rref)

V = cat(1, lref, rref)';

%coords = rc * inv(V);
coords = rc / V;

if abs(coords(:,1)) == -coords(:,1)
    coords(:,1) = -coords(:,1);
    l_flip = true;
else
    l_flip = false;
end

if abs(coords(:,2)) == -coords(:,2)
    coords(:,2) = -coords(:,2);
    r_flip = true;
else
    r_flip = false;
end

[junk i_corner] = min(sum(coords, 2));  %#ok<ASGLU>

end

% 
% function squareFactors = makeSquare(n)
% 
% rotMat = [0 1; -1 0];
% 
% v = zeros(2 * n, 2);
% v(:,2) = n;
% v(:,1) = -(n - 1):n;
% 
% squareFactors = v;
% 
% for ii = 1:3
%     v = v * rotMat;
%     squareFactors = cat(1, squareFactors, v);
% end
% end