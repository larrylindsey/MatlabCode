function [R, T] = least_squares_rigidified(match_points_a, match_points_b, weights, ...
  fixed_patch_id) %#ok<INUSL>
% [R, T] = least_squares_rigidified(match_points_a, match_points_b, weights, ...
%   fixed_patch_id)
% does a least-squares fit of a file of transform points. Rigidifies the
% transforms to have scale nearly 1.
%
% Input:
%   match_points_a    NxN cell matrix where N is the number of patches.
%   match_points_b    NxN cell matrix where N is the number of patches.
%     The ith pair of correspondence points from patch p1 to p2 is given by
%     x: match_points_a{p1,p2}(1,i) <--> match_points_b{p1,p2}(1,i)
%     y: match_points_a{p1,p2}(2,i) <--> match_points_b{p1,p2}(2,i)
%   weights           NxN double matrix of weights to be given for
%                       corresponding patch p1 to p2.
%   fixed_patch_id    patch that should be given identity transformation.
%
% Output:
%   R         1xN cell array of rotation-skew part of the affine transform
%               A(1:2,1:2). R{p} corresponds to patch p.
%   T         1xN cell array of translation part of the affine transform
%               A(1:2,3). T{p} corresponds to patch p.
%
% Lou Scheffer
% Fellow,
% Janelia Farm Research Campus, HHMI.
% July 1st, 2009.
%

global alignLHS alignRHS
% [n1,x1,y1,n2,x2,y2] = textread(name,'%d %f %f %d %f %f');
% if max(n1) > max(n2)
%   nt = max(n1)+1;
% else
%   nt = max(n2)+1;
% end

% Now nt= number of transforms we need to find
n_patch = size(match_points_a,1);
scale = 10000;  % scale coordinates down by this
nvars = n_patch*6; %each transform has 6 variables we need to set

alignLHS = sparse([],[],[],nvars,nvars,100); alignRHS = zeros(nvars, 1);  %% allow 100 non-zero terms to start

fprintf('Setting up the equations ...\n');
tic
for patch_1 = 1:n_patch
  t1 = patch_1 - 1; % the transform id - 0 based
  j = t1*6; % base variable number
  for patch_2 = 1:n_patch
    if(isempty(match_points_a{patch_1, patch_2}))
      continue;
    end
    fprintf('patch_1: %d, patch_2: %d\n', patch_1, patch_2);
    
    t2 = patch_2 - 1; % the transform id - 0 based
    k = t2 * 6; % base variable number
    
    %%% Add in the x-coords correspondences
    % the indexes of the variables within the constraint matrix used in
    % x-coord least squares
    indices = [j+1,j+2,j+3,k+1,k+2,k+3];
    for match_id = 1:size(match_points_a{patch_1, patch_2}, 2)
      vals = [match_points_a{patch_1, patch_2}(1,match_id), ...
        match_points_a{patch_1, patch_2}(2,match_id), scale, ...
        -match_points_b{patch_1, patch_2}(1,match_id), ...
        -match_points_b{patch_1, patch_2}(2,match_id), -scale]/scale;
      AddConstraint(indices, vals, 0.0);
    end

    %%% Add in the y-coords correspondences
    % the indexes of the variables within the constraint matrix used in
    % y-coord least squares
    indices = [j+4,j+5,j+6,k+4,k+5,k+6];
    for match_id = 1:size(match_points_a{patch_1, patch_2}, 2)
      vals = [match_points_a{patch_1, patch_2}(1,match_id), ...
        match_points_a{patch_1, patch_2}(2,match_id), scale, ...
        -match_points_b{patch_1, patch_2}(1,match_id), ...
        -match_points_b{patch_1, patch_2}(2,match_id), -scale]/scale;
      AddConstraint(indices, vals, 0.0);
    end
  end
end

fprintf('Setting patch %d with identity transform ...\n', fixed_patch_id);
%Now add the constraints for the initial square
InitStiff = 10000/scale;
j = (fixed_patch_id-1)*6;
AddConstraint(j+1, InitStiff, InitStiff);
AddConstraint(j+2, InitStiff, 0.0);
AddConstraint(j+3, InitStiff, 0.0);
AddConstraint(j+4, InitStiff, 0.0);
AddConstraint(j+5, InitStiff, InitStiff);
AddConstraint(j+6, InitStiff, 0.0);

%controls stiffness by which we seek similarity
fprintf('Setting stiffness constraints for the first pass ...\n');
stiff_sim = 10000/scale;  
for i = 1:n_patch
  j=(i-1)*6;   %base variable number
  AddConstraint([j+1,j+5], [stiff_sim, -stiff_sim], 0.0);
  AddConstraint([j+2,j+4], [stiff_sim,  stiff_sim], 0.0);
  
end
toc

fprintf('--Starting first solve...\n');
tic
X = alignLHS\alignRHS;
if(fixed_patch_id>n_patch/2)
  probe_patch = 1;
else
  probe_patch = n_patch;
end
toc
fprintf('Transform of patch farthest from identity patch:\n');
fprintf('%f ',X(probe_patch*6-5:probe_patch*6)); fprintf('\n');

fprintf('Adding constraints to control scale (circle tangent) ...\n;');
tic
scale_stiff = 10000/scale;
for i=1:n_patch
  j = (i-1)*6;  % j is the base of the transform's 6 variables
  a = X(j+1);
  b = X(j+2);
  m = sqrt(a*a+b*b);
  %a,b,m
  AddConstraint([j+1, j+2], [a/m*scale_stiff, b/m*scale_stiff], scale_stiff);
end
toc

fprintf('--Starting second solve...\n');
tic
X = alignLHS\alignRHS;
toc

% Scale the translation coordinates back to original resolution
for i=1:n_patch
  j = (i-1)*6;
  X(j+3) = X(j+3)*scale;
  X(j+6) = X(j+6)*scale;
end
fprintf('Transform of patch farthest from identity patch:\n');
fprintf('%f ',X(probe_patch*6-5:probe_patch*6)); fprintf('\n');

%calculate the residuals
fprintf('Calculating residuals and putting the transforms for output ...\n');
total_residual_sqr = 0.0;
n_total_match = 0;
R = cell(1, n_patch);
T = cell(1, n_patch);
for patch_1 = 1:n_patch
  t1 = patch_1 - 1; % the transform id - 0 based
  j = t1*6; % base variable number
  
  for patch_2 = 1:n_patch
    if(isempty(match_points_a{patch_1, patch_2}))
      continue;
    end
    fprintf('patch_1: %d, patch_2: %d\n', patch_1, patch_2);
    
    t2 = patch_2 - 1; % the transform id - 0 based
    k = t2 * 6; % base variable number
    
    A = [X(j+1), X(j+2), X(j+3); X(j+4), X(j+5), X(j+6)];
    p_t_1 = A*[match_points_a{patch_1, patch_2}; ...
      ones(1, size(match_points_a{patch_1, patch_2},2))];
    
    A = [X(k+1), X(k+2), X(k+3); X(k+4), X(k+5), X(k+6)];
    p_t_2 = A*[match_points_b{patch_1, patch_2}; ...
      ones(1, size(match_points_b{patch_1, patch_2},2))];
    
    error_p = p_t_1 - p_t_2;
    total_residual_sqr = sum(sum(error_p.^2));
    n_total_match = n_total_match + size(error_p, 2);
    
  end
  
  R{patch_1} = [X(j+1), X(j+2); X(j+4), X(j+5)];
  T{patch_1} = [X(j+3); X(j+6)];
end
residual_rms = sqrt(total_residual_sqr / n_total_match);
fprintf('RMS error %f pixels\n', residual_rms);

return
end

