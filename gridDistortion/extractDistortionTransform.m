function [rc_found, rc_match_aff, rc_grid, dbStr] =...
    extractDistortionTransform(im0, varargin)
tic;
if ischar(im0)
    im0 = imread(im0);
end

if size(im0, 3) > 1
    im0 = rgb2gray(im0);
end

if ~isfloat(im0)
    im0 = im2single(im0);
end


HH = [];
match = [];
rStr = [];
bwEnergy = [];
scale = size(im0);
scale = reshape(scale(1:2), [1 2]);

while ~isempty(varargin)
    arg = varargin{1};
    i_match = strmatch(arg, {'HH', 'match', 'rStr', 'bwEnergy'});%#ok
    i_match = i_match(1);
    switch i_match
        case 1
            disp('Got HH');
            HH = varargin{2};
        case 2
            disp('Got match');
            match = varargin{2};
        case 3
            disp('Got rStr');
            rStr = varargin{2};
        case 4
            disp('Got bwEnergy');
            bwEnergy = varargin{2};
        otherwise
            error('Unexpected value');
    end
    varargin = varargin(3:end);
end

%Step 0: Fix orientation of image
s = cputime;
if isempty(rStr) && isempty(match)
    fprintf('Finding Approximate Grid\n');
    rStr = findRoughLines(im0);
    save -v7.3 grid_approx_cache.mat rStr    
end
s(end + 1) = cputime;

if isempty(match)
    fprintf('Creating match filter\n');

    match = gridEstimate(rStr, im0);
    if nargout > 3
        dbStr.rStr = rStr;
    end
    cm = rStr.cropMask;
    clear rStr;
    rStr.cropMask = cm;
    clear cm;
    save match_cache.mat match

    imwrite(match, 'match_attach.png');
end
s(end + 1) = cputime;

%Step 1: Label Grid Intersections.
sprintf('Labeling Grid Intersections\n');
if isempty(HH)

    %[rmin rmax cmin cmax] = mask2boundingbox(rStr.cropMask);
    %im0(not(rStr.cropMask)) = 0;
%     im0 = im0(rmin:rmax, cmin:cmax);
    [rc_found, grid_model, model_err, L, HH] =...
        labelGridIntersections(im0, match, []);
%     rc_found(:,1) = rc_found(:,1) + rmin - 1;
%     rc_found(:,2) = rc_found(:,2) + cmin - 1;

    save -v7.3 match_energy_cache.mat HH

elseif isempty(bwEnergy)

    [rc_found, grid_model, model_err, L, HH] =...
        labelGridIntersections(im0, match, [], HH);

else

    [rc_found, grid_model, model_err, L, HH] =...
        labelGridIntersections(im0, match, [], HH, bwEnergy);
    
end
s(end + 1) = cputime;

if nargout > 1
    dbstop if error;
    dbStr.HH = HH;
end

clear HH;

%[rc_grid grid_model] = fitGrid(grid_model, rc_found, model_err);

rc_grid = matchGridToSquareGrid(rc_found, grid_model, ...
    sqrt(abs(det(grid_model))) / 10);
rc_found = rc_found(rc_grid(:,1), :);
rc_grid = rc_grid(:,[2 3]);

rc_found = rc_found ./ repmat(scale / 2, [size(rc_found,1), 1]) - 1;
grid_model = grid_model ./ repmat(scale / 2, [2 1]);

[rc_grid grid_model] = fixGridMatch(rc_grid, grid_model);
rc_grid = rc_grid(:,[2 1]);

%grid_model = eye(2) * sqrt(abs(det(grid_model)));

% rc_match_sq = getSquareMatchGrid(rc_found, getSquareGM(grid_model), rc_grid);
rc_match_aff = getAffineMatchGrid(rc_found, getSquareGM(grid_model), rc_grid);


if nargout > 1
    dbStr.L = L;
    dbStr.match = match;
    dbStr.s = s;
    dbStr.rc_found = rc_found;
    dbStr.rc_grid = rc_grid;
    dbStr.grid_model = grid_model;
    dbStr.model_err = model_err;
end

% c = clock;
% sendmsg(msgsubject,...
%     sprintf('Finished processing at %d:%d:%g\n', c(4), c(5), c(6)));
toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rc_grid grid_model] = fixGridMatch(rc_grid, grid_model)

e1_ind = abs(grid_model(:,1));
if e1_ind(2) > e1_ind(1)
    rc_grid = rc_grid(:,[2 1]);
    grid_model = grid_model(:,[2 1]);
end

sign_ind = grid_model(logical(eye(2)));

if sign_ind(1) < 0
    grid_model(1,:) = - grid_model(1,:);
    rc_grid(:,1) = -rc_grid(:,1);
end

if sign_ind(2) < 0
    grid_model(2,:) = -grid_model(2,:);
    rc_grid(:,2) = -rc_grid(:,2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gmSQ = getSquareGM(grid_model)

gmDet = det(grid_model);
d = sqrt(abs(gmDet));

gmSQ = eye(2) * d;

if gmDet < 0
    warning('extractedTransform:swappedAxes',...
        'Detected swapped axes in grid model');
    gmSQ = gmSQ([2 1],:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % function matchGrid = getSquareMatchGrid(rc_pt, gm, rc_grid)
% % 
% % grid0 = rc_grid * gm;
% % matchGrid = trsAlign(grid0, rc_pt);
% % 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function matchGrid = getAffineMatchGrid(rc_pt, gm, rc_grid)

grid0 = rc_grid * gm;
tr = fitTransform(grid0, rc_pt, 1);
matchGrid = applyTransform(grid0, tr);


end
