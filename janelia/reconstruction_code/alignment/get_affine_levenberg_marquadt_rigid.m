function [a,b] = get_affine_levenberg_marquadt_rigid(...
  match_points_a, match_points_b, constant_patch_id, lm_config)
% [a,b] = get_affine_levenberg_marquadt_rigid(match_points_1_a, ...
%           match_points_1_b)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  01142009  init. code
%

% Notation:
% beta: vector of affine parameters for all the patches
% y: different from the coordinates. Here y refers to generic desired
% function value

n_patch = size(match_points_a, 1);

% scaling parameters in LM
nu = 2; % lambda scaling factor
lambda = lm_config.lambda; % 0.0001;
alpha_1 = lm_config.alpha; % 1000;
alpha_2 = lm_config.alpha; % 1000;
alpha_3 = lm_config.alpha; % 1000;
gamma = lm_config.gamma; % 0.1;
max_n_iteration = lm_config.max_n_iteration; % 1000;

% the patch to kept constant with identity transformation
forward_ref = [1:(constant_patch_id-1), 0, (constant_patch_id):(n_patch-1)];

% initial guess
% % Use linear least-squares.
% fprintf('Computing initial guess using least-squares ... ');
% [a_init,b_init] = getSIFTmodel(match_points_a, match_points_b, ...
%   ones(size(match_points_a)), [], constant_patch_id);
% beta = zeros(6 * (n_patch-1), 1);
% i=0;
% for j = 1:length(a_init)
%   if(j==constant_patch_id)
%     continue;
%   end
%   %   a{patch_id} = beta_s([1 2;4 5]);
%   %   b{patch_id} = beta_s([3 6]);
%   beta(i*6+1:i*6+6) = [a_init{j}(1,1), a_init{j}(1,2), b_init{j}(1), ...
%     a_init{j}(2,1), a_init{j}(2,2), b_init{j}(2)];
%   i = i+1;
% end
% fprintf('done.\n');
% % Use random
beta = 0.5+0.25*rand(6 * (n_patch-1), 1);

% randomization mask
random_mask_1 = zeros(size(beta));
random_mask_1(3:6:end) = 1;
random_mask_1(6:6:end) = 1;

random_mask_2 = zeros(size(beta));
random_mask_2(1:6:end) = 1;
random_mask_2(2:6:end) = 1;
random_mask_2(4:6:end) = 1;
random_mask_2(5:6:end) = 1;

% Levenberg-Marquadt iterations
delta_vals = [];
fprintf('Precomputing constant part of the Jacobian ...\n');
n_row = 0;
for patch_id_1 = 1:n_patch
  for patch_id_2 = patch_id_1+1:n_patch
    if(isempty(match_points_a{patch_id_1, patch_id_2}))
      continue;
    end
    n_row = n_row + 2*size(match_points_b{patch_id_1, patch_id_2},2);
  end
end
J_constant = zeros(n_row, 6*(n_patch-1));
row_id = 0;
for patch_id_1 = 1:n_patch
  fprintf('%d ', patch_id_1);
  if(mod(patch_id_1, 20)==0)
    fprintf('\n');
  end
  for patch_id_2 = patch_id_1+1:n_patch
    if(isempty(match_points_a{patch_id_1, patch_id_2}))
      continue;
    end
    J_temp = zeros(2*size(match_points_b{patch_id_1, patch_id_2},2), 6*(n_patch-1));
    if(patch_id_1~=constant_patch_id)
      J_temp(1:2:end, ...
        (forward_ref(patch_id_1)-1)*6+1:(forward_ref(patch_id_1)-1)*6+3) = ...
        [match_points_a{patch_id_1, patch_id_2}', ...
        ones(size(match_points_a{patch_id_1, patch_id_2},2),1)];
      J_temp(2:2:end, ...
        (forward_ref(patch_id_1)-1)*6+4:forward_ref(patch_id_1)*6) = ...
        [match_points_a{patch_id_1, patch_id_2}', ...
        ones(size(match_points_a{patch_id_1, patch_id_2},2),1)];
    end
    if(patch_id_2~=constant_patch_id)
      J_temp(1:2:end, ...
        (forward_ref(patch_id_2)-1)*6+1:(forward_ref(patch_id_2)-1)*6+3) = ...
         -  [match_points_b{patch_id_1, patch_id_2}', ...
        ones(size(match_points_b{patch_id_1, patch_id_2},2),1)];
      J_temp(2:2:end, ...
        (forward_ref(patch_id_2)-1)*6+4:forward_ref(patch_id_2)*6) = ...
         -  [match_points_b{patch_id_1, patch_id_2}', ...
        ones(size(match_points_b{patch_id_1, patch_id_2},2),1)];
    end

    J_constant(row_id+1:row_id+size(J_temp,1),:) = J_temp;
    row_id = row_id + size(J_temp,1);
  end
end

% Compute sizes of datastructure for pre-allocation
size_delta_y_corresp = [0, 1];
for patch_id_1 = 1:n_patch
  beta_1 = get_beta(patch_id_1);
  tform_1 = beta_1([1 2 3; 4 5 6]);
  for patch_id_2 = patch_id_1+1:n_patch
    if(isempty(match_points_a{patch_id_1, patch_id_2}))
      continue;
    end
    size_delta_y_corresp(1) = size_delta_y_corresp(1) + ...
      2*size(match_points_a{patch_id_1, patch_id_2},2);
  end
end

size_J_regularizer = [0, 6*(n_patch-1)];
n_skip_row = 1;
for patch_id = 1:n_skip_row:n_patch
  if(patch_id==constant_patch_id)
    continue;
  end
  size_J_regularizer(1) = size_J_regularizer(1) + 3;
end

fprintf('for the correspondences\n');
tic
row_id = [];
beta_1 = [];
beta_2 = [];
tform_1 = [];
tform_2 = [];
ones_temp = [];
delta_y_corresp = get_delta_y_corresp();
toc

J_regularizer = zeros(size_J_regularizer);
delta_y_regularizer = zeros(size_J_regularizer(1),1);
for iter = 1:max_n_iteration
  fprintf('for the regularization\n');
  tic
  row_id = 1;
  for patch_id = 1:n_skip_row:n_patch
    if(patch_id==constant_patch_id)
      continue;
    end
    beta_s = get_beta(patch_id);
    % ab + cd = 0
    delta_y(row_id) = ...
      alpha_1 * (- (beta_s(1)*beta_s(2) + beta_s(4)*beta_s(5)));
    J_regularizer(row_id, ...
      (forward_ref(patch_id)-1)*6+1:forward_ref(patch_id)*6) = ...
      alpha_1 * [beta_s(2), beta_s(1), 0, beta_s(5), beta_s(4), 0];
    % ad - bc = 1
    delta_y(row_id+1) = ...
      alpha_2 * (1 - (beta_s(1)*beta_s(5) - beta_s(2)*beta_s(4)));
    J_regularizer(row_id+1, ...
      (forward_ref(patch_id)-1)*6+1:forward_ref(patch_id)*6) = ...
      alpha_2 * [beta_s(5), -beta_s(4), 0, -beta_s(2), beta_s(1), 0];
    % a^2 + c^2 - b^2 - d^2 = 0
    delta_y(row_id+2) = ...
      alpha_3 * (- (beta_s(1)^2 + beta_s(4)^2 - beta_s(2)^2 - beta_s(5)^2));
    J_regularizer(row_id+2, ...
      (forward_ref(patch_id)-1)*6+1:forward_ref(patch_id)*6) = ...
      alpha_3 * [2*beta_s(1), -2*beta_s(2), 0, 2*beta_s(4), -2*beta_s(5), 0];
    row_id = row_id + 3;
  end
  if(nnz(size(delta_y_corresp)~=size_delta_y_corresp)~=0)
    error('change in size');
  end
  if(nnz(size(J_regularizer)~=size_J_regularizer)~=0)
    error('change in size');
  end
  toc
  delta_y = [delta_y_corresp; delta_y_regularizer];
  delta_y = reshape(delta_y, [length(delta_y), 1]);
  
  J = [J_constant; J_regularizer];

  H = J'*J;
  R = J'*delta_y;

  L = H + lambda*diag(diag(H));
  fprintf('linsolve 1\n');
  tic
  delta_1 = linsolve(L, R);
  toc

  L = H + (lambda/nu)*diag(diag(H));
  fprintf('linsolve 2\n');
  tic
  delta_2 = linsolve(L, R);
  toc

  beta_old = beta;

  beta = beta_old + gamma*delta_1;
  fprintf('for the correspondences\n');
  tic
  delta_y_corresp = get_delta_y_corresp();
  toc
  delta_y_1 = delta_y_corresp;
  error_1 = sqrt(mean(delta_y_1.^2));

  beta = beta_old + gamma*delta_2;
  fprintf('for the correspondences\n');
  tic
  delta_y_corresp = get_delta_y_corresp();
  toc
  delta_y_2 = delta_y_corresp;
  error_2 = sqrt(mean(delta_y_2.^2));

  if(error_2<error_1)
    lambda = lambda/nu;
    err = error_2;
    beta = beta_old + gamma*delta_2;
    delta_y_corresp = delta_y_2;
  else
    if(~isempty(delta_vals) && error_1>delta_vals(end))
      lambda = lambda*nu;
      beta_rand = 2*random_mask_1.*(rand(size(random_mask_1))-0.5) + ...
          0.1*random_mask_2.*(rand(size(random_mask_2))-0.5);
      beta = beta_old;% + beta_rand;
      err = error_1;
      delta_y_corresp = delta_y_1;
    else
      beta = beta + gamma*delta_1;
      err = error_1;
      delta_y_corresp = delta_y_1;
    end
  end

  delta_vals(end+1) = err;

  if(lm_config.is_verbose)
    fprintf('iteration: %d, max err: %d, RMS err: %d\n', iter, ...
            max(abs(delta_y_corresp)), delta_vals(end));
  end
  if(length(delta_vals)>1)
    if(delta_vals(end-1) > delta_vals(end) && ...
          (delta_vals(end-1) - delta_vals(end)) < -lm_config.min_error_change)
      break;
    end
  end
end

if(lm_config.is_verbose_figures)
  figure(1); plot(log(delta_vals));
end

a = cell(n_patch, 1);
b = cell(n_patch, 1);
for patch_id = 1:n_patch
  beta_s = get_beta(patch_id);
  a{patch_id} = beta_s([1 2;4 5]);
  b{patch_id} = beta_s([3 6]);
end
return

  function beta_s = get_beta(patch_id)
    if(patch_id == constant_patch_id)
      beta_s = [1 0 0 0 1 0]';
    else
      beta_s = beta((forward_ref(patch_id)-1)*6+1:forward_ref(patch_id)*6);
    end
    return
  end

  function delta_y_corresp = get_delta_y_corresp()
    delta_y_corresp = zeros(size_delta_y_corresp); % error
    row_id = 0;
    for patch_id_1 = 1:n_patch
      beta_1 = get_beta(patch_id_1);
      tform_1 = beta_1([1 2 3; 4 5 6]);
      for patch_id_2 = patch_id_1+1:n_patch
        if(isempty(match_points_a{patch_id_1, patch_id_2}))
          continue;
        end
        beta_2 = get_beta(patch_id_2);
        tform_2 = beta_2([1 2 3; 4 5 6]);
        ones_temp = ones(1, size(match_points_a{patch_id_1, patch_id_2},2));
        delta_y_corresp(...
          row_id+1:row_id+2*size(match_points_a{patch_id_1, patch_id_2},2))...
          = reshape(...
          - tform_1 * [match_points_a{patch_id_1, patch_id_2}; ones_temp] + ...
          tform_2 * [match_points_b{patch_id_1, patch_id_2}; ones_temp], ...
          [1, 2*size(match_points_a{patch_id_1, patch_id_2},2)]);
        row_id = row_id + 2*size(match_points_a{patch_id_1, patch_id_2},2);
      end
    end
    return
  end
end
