function [rc_found, rc_grid, grid_model, dbStr] =...
    extractDistortionTransform(im0, varargin)

if ischar(im0)
    im0 = imread(im0);
end

if size(im0, 3) > 1
    im0 = rgb2gray(im0);
end
im0 = im2single(im0);


HH = [];
% imGrid = [];
match = [];
rStr = [];
bwEnergy = [];

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
    varargin = {varargin{3:end}};
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
    [rmask cmask] = find(rStr.cropMask);
    rmin = min(rmask); rmax = max(rmask);
    cmin = min(cmask); cmax = max(cmask);
    im0 = im0 .* rStr.cropMask;
    im0 = im0(rmin:rmax, cmin:cmax);
    [rc_found, grid_model, model_err, L, HH] =...
        labelGridIntersections(im0, match,...
        rStr.cropMask(rmin:rmax,cmin:cmax));
    rc_found(:,1) = rc_found(:,1) + rmin - 1;
    rc_found(:,2) = rc_found(:,2) + cmin - 1;
    sendmsg(msgsubject, 'Labeled Grid Intersections.');
    save match_energy_cache.mat HH
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

[rc_grid grid_model] = fitGrid(grid_model, rc_found, model_err);


% fprintf('Calculating Match Grid\n');
% tic;
% matchGrid = fftfilter2(imGrid, match);
% toc;
% 
% sendmsg(msgsubject, 'Calculated Match Grid');
% 
% s(end + 1) = cputime;
% 
% %Step 2: Use the grid intersection labels to find the exact grid
% %intersection points.  This is noisy.
% fprintf('Finding cross points\n');
% tic;
% [r c] = findCrossPoints(matchGrid, L);
% [rlReg clReg] = findBetterRegularCrossPoints(r, c, rowModel, colModel);
% %[rlReg clReg] = findRegularCrossPoints(r, c, rl, cl);
% toc;
% s(end + 1) = cputime;
% 
% 
% %Step 3: Calculate the deviation from a linear grid.
% fprintf('Calculating grid nonlinearity\n');
% tic;
% [rTest cTest] = meshgrid(rlReg, clReg);
% %dv = gridDifference(r, c, rTest(:), cTest(:));
% [rTest cTest] = assignGrid(r, c, rTest, cTest);
% % Now rTest(p) is the closest "rectangular grid" point to r(p), for a given
% % index p.
% sendmsg(msgsubject, 'Found cross points and grid nonlinearity');
% 
% toc;
% s(end + 1) = cputime;
% 
% %Step 4: Calculate the transform that will remove the 2nd order
% %nonlinearities.
% fprintf('Calculating Transform\n');
% tic;
% [fFwd, fInv, T, transinfo] = getTransformInfo(r, c, rTest, cTest);
% toc;
% s(end + 1) = cputime;
% 
% %Now create the transform object that will be used with imtform.
% fprintf('Creating transform struct\n');
% tic;
% Tr = maketform('custom', 2, 2, fFwd, fInv, T);
% s(end + 1) = cputime;
% toc;
% fprintf('Total computation time: %f\n', s(end) - s(1));

if nargout > 1
    dbStr.L = L;
    dbStr.match = match;
    dbStr.s = s;
    dbStr.rc_found = rc_found;
    dbStr.rc_grid = rc_grid;
    dbStr.grid_model = grid_model;
end

c = clock;
sendmsg(msgsubject,...
    sprintf('Finished processing at %d:%d:%g\n', c(4), c(5), c(6)));

end