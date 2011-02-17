function startup()
% boot up the EM reconstruction pipeline
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

clear all
close all

fprintf('--------------------------------------------------\n');
fprintf('| EM Reconstruction and Visualization Software   |\n');
fprintf('| Chklovskii Lab., Janelia Farm Research Campus, |\n');
fprintf('| Howard Hughes Medical Institute                |\n');
fprintf('--------------------------------------------------\n\n');

%%%
% Set MATLAB paths
%%%
global code_dir

fprintf('Setting MATLAB paths for EM reconstruction pipeline ...');
scripts_dir = [strrep(pwd(), '\', '/'), '/'];

% check to make sure we are in the correct dir.
dir_brackets = strfind(scripts_dir, '/');
if(strcmp(scripts_dir(dir_brackets(end-2):end), '/reconstructions/scripts/')~=1)
  error(['Not in correct directory. Please start MATLAB in ', ...
    'em_reconstruction/reconstructions/scripts/']);
end

% check if developer or user
is_developer = false;
dev_code_dir = [scripts_dir(1:dir_brackets(end-2)), 'code/'];
if(exist(dev_code_dir, 'dir')==7)
  is_developer = true;
end

% set code directory
if(is_developer)
  code_dir = dev_code_dir;
else
  code_dir = '/groups/chklovskii/chklovskiilab/em_reconstruction/code/';
end
addpath(code_dir);

% set paths
set_matlab_path(code_dir);

% set some warnings off
warning('OFF', 'Images:initSize:adjustingMag');

fprintf('done.\n');
return
end
