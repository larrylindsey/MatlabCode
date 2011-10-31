function [rc_found, rc_grid, grid_model, scale, dbStr] =...
    extractDistortionTransform(im0, varargin)

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
% imGrid = [];
match = [];
rStr = [];
bwEnergy = [];
scale = size(im0);
scale = reshape(scale(1:2), [1 2]);

while ~isempty(varargin)
    arg = varargin{1};
    i_match = strmatch(arg, {'HH', 'imGrid', 'match', 'rStr', 'bwEnergy'});
    i_match = i_match(1);
    switch i_match
        case 1
            disp('Got HH');
            HH = varargin{2};
        case 2
            disp('Got imGrid');
%             imGrid = varargin{2};
        case 3
            disp('Got match');
            match = varargin{2};
        case 4
            disp('Got rStr');
            rStr = varargin{2};
        case 5
            disp('Got bwEnergy');
            bwEnergy = varargin{2};
        otherwise
            error('Unexpected value');
    end
    varargin = varargin(3:end);
end


msgsubjectbeg = 'Distortion Extraction Output';
msgsubject = ['RE: ' msgsubjectbeg];
%Step 0: Fix orientation of image
s = cputime;
if isempty(rStr)
    tic;
    fprintf('Finding Approximate Grid\n');
    rStr = findRoughLines(im0);
    toc;
    sendmsg(msgsubjectbeg, 'Found Rough Grid Approximation');
    save -v7.3 grid_approx_cache.mat rStr    
end

if isempty(match)
    fprintf('Creating match filter\n');
    tic;
    match = gridEstimate(rStr, im0);
    cm = rStr.cropMask;
    clear rStr;
    rStr.cropMask = cm;
    clear cm;
    save match_cache.mat match
    toc;
    sendmsg(msgsubjectbeg, ...
        sprintf('Created Match Filter of size %d by %d', size(match, 1), ...
        size(match, 2)));
    imwrite(match, 'match_attach.png');
    try
        sendmail('larry.f.lindsey@gmail.com', 'Match Kernel', 'Attached', ...
            'match_attach.png');
    catch senderr
        fprintf('Unable to send match email.  Error was %s\n',...
            senderr.message);
    end
end

% save -v7.3 grid_stuff rStr match imGrid
% 
% keyboard;

%Step 1: Label Grid Intersections.
sprintf('Labeling Grid Intersections\n');
if isempty(HH)
    tic;
    [rmin rmax cmin cmax] = mask2boundingbox(rStr.cropMask);
%     [rmask cmask] = find(rStr.cropMask);
%     rmin = min(rmask); rmax = max(rmask);
%     cmin = min(cmask); cmax = max(cmask);
    im0(not(rStr.cropMask)) = 0;
    %im0 = im0 .* rStr.cropMask;
    im0 = im0(rmin:rmax, cmin:cmax);
    [rc_found, grid_model, model_err, L, HH] =...
        labelGridIntersections(im0, match,...
        rStr.cropMask(rmin:rmax,cmin:cmax));
    rc_found(:,1) = rc_found(:,1) + rmin - 1;
    rc_found(:,2) = rc_found(:,2) + cmin - 1;
    sendmsg(msgsubject, 'Labeled Grid Intersections.');
    save -v7.3 match_energy_cache.mat HH
%    keyboard;
%    close all;
%    imagesc(HH);
%    drawnow;
    toc;
elseif isempty(bwEnergy)
    tic;
    [rc_found, grid_model, model_err, L, HH] =...
        labelGridIntersections(im0, match, rStr.cropMask, HH);
    toc;
else
    tic;
    [rc_found, grid_model, model_err, L, HH] =...
        labelGridIntersections(im0, match, rStr.cropMask, HH, bwEnergy);
    toc;
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

if nargout > 1
    dbStr.L = L;
    dbStr.match = match;
    dbStr.s = s;
    dbStr.rc_found = rc_found;
    dbStr.rc_grid = rc_grid;
    dbStr.grid_model = grid_model;
    dbStr.model_err = model_err;
end

c = clock;
sendmsg(msgsubject,...
    sprintf('Finished processing at %d:%d:%g\n', c(4), c(5), c(6)));

end

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