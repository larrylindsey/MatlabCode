function start_proofread_session2(config_in)
% Entry point into proofreading. Execute this to start any proofreading session.
%
% When proofreading first go to <chklovskii_lab>/em_reconstruction
% directory and then launch Matlab. This would ensure correct path
% settings.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04082008  init code
% v1  04082008  GUI dialog boxes
% v2  04142008  inclusion of reconstruction pipeline arch.
% v3  04162008  changing the config structure; to be launched from
%               proofreading directory.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear al cat superpixel_2_seg_map link_threshold proof data debug
close all

if(nargin<1 || isempty(config_in))
  config.to_be_proofread.proofreading_sessions_dir = ...
    '/groups/chklovskii/chklovskiilab/em_reconstruction/data_to_be_proofread/';
  config.proofreading_session.proofreading_sessions_dir = ...
    '/groups/chklovskii/chklovskiilab/em_reconstruction/proofreading_sessions/';
else
  config.to_be_proofread.proofreading_sessions_dir = ...
    config_in.to_be_proofread.root_dir;
  config.proofreading_session.proofreading_sessions_dir = ...
    config_in.proofreading_session.root_dir;
end

tobe_proofread_dir = config.to_be_proofread.proofreading_sessions_dir;
proofreading_sessions_dir = config.proofreading_session.proofreading_sessions_dir;

[status, user_name] = system('id -u -n');
if(status==0)
  user_name = user_name(1:end-1);
else
  user_name = input('Enter your user-name: ', 's');
end

fprintf('-------------------------------------\n');
fprintf(['Hello ', user_name, '\n']);
fprintf('Starting a proofreading session\n');
fprintf('-------------------------------------\n\n');

if(exist('_check_path_em_reconstruction_code.m', 'file')~=2)
  error('The path to code directories not set, exiting');
end;

cd(proofreading_sessions_dir);
if(exist(user_name, 'dir')~=7)
  mkdir2(user_name);
end;
proofreading_sessions_dir = [proofreading_sessions_dir, user_name, '/'];
cd(user_name);

global al cat superpixel_2_seg_map links_3D reconstruction_config

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step I: Get the proofreading session information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist('./_start_proofread_session.log.mat', 'file')==2) % give previous session as an option for quick launch
  prev_session_info = load2('./_start_proofread_session.log.mat');
end;

options = cell(3,1);
options{1} = 'Continue an existing session';
options{2} = 'Copy an existing session to a new session and start from there';
options{3} = 'Start a new proofreading session';
if(exist('prev_session_info', 'var')==1)
  options{4} = ['Continue previous session - name: ', prev_session_info.name];
end;
[o, v] = listdlg('PromptString', 'Choose an option', 'SelectionMode', 'single', 'ListString', options, 'ListSize', [500 80]);
if(ismember(v, 1:length(options))==0)
  error('Could not determine option, exiting');
end;

switch(o)
  case 1 % continue an existing session
    cd(proofreading_sessions_dir);
    [session_file_name, session_path_name] = uigetfile('*.session_info.mat', 'Choose a session file');
    cd(proofreading_sessions_dir);
    if(session_file_name==0)
      error('Could not determine session, exiting.\n');
    end;
    
    % load the session information
    session_info = load2([session_path_name, session_file_name]);
    
  case 2 % Copy an existing session to a new session and start from there
    cd(proofreading_sessions_dir);
    [existing_session_file_name, existing_session_path_name] = uigetfile('*.session_info.mat', 'Choose a session file');
    cd(proofreading_sessions_dir);
    if(existing_session_file_name==0)
      error('Could not determine session, exiting.\n');
    end;
    existing_session_info = load2([existing_session_path_name, existing_session_file_name]);

    % get user information regarding the new copy of this session
    session_info = existing_session_info;
    op.WindowStyle = 'normal';
    session_info.name = char(inputdlg('Naming convention <proofreader name>.<date>.<additional_info>.', ...
      'Give session name', 1, {''}, op));

    % make a new directory for the new session
    session_path_name = [proofreading_sessions_dir, session_info.stack_name, '/', ...
      session_info.region_name, '/', session_info.reconstruction_name, '/', ...
      session_info.name, '/'];
    mkdir2(session_path_name);
    system(['chmod a+rw ', session_path_name]);
    
    % copy files to the new session
    fprintf('\nCopying existing session data to new session ..');
%     system(['cp -r ', existing_session_path_name, '/* ', session_path_name, '/']);
    copyfile([existing_session_path_name, '/*'], [session_path_name, '/']);
%     system(['rm -rf ', session_path_name, '/*.session_info.mat']);
    delete([session_path_name, '/*.session_info.mat']);

    % save for later access
    save2([session_path_name, session_info.name, '.session_info.mat'], '-STRUCT', 'session_info');

    fprintf('done\n');
    
  case 3 % Start a new session
    session_info = [];
    cd(tobe_proofread_dir);
    curr_dir = pwd2;
    stack_dir = uigetdir('.', 'Give the stack name');
    session_info.stack_name = stack_dir(length(curr_dir)+2:end);

    cd(stack_dir);
    curr_dir = pwd2;
    region_dir = uigetdir('.', 'Give the region name');
    session_info.region_name = region_dir(length(curr_dir)+2:end);

    cd(region_dir);
    curr_dir = pwd2;
    reconstruction_dir = uigetdir('.', 'Give the reconstruction name');
    session_info.reconstruction_name = reconstruction_dir(length(curr_dir)+2:end);
    cd(proofreading_sessions_dir);
    
    op.WindowStyle = 'normal';
    session_info.name = char(inputdlg('Naming convention <proofreader name>.<date>.<additional_info>.', ...
      'Give session name', 1, {''}, op));
    
    fprintf('Enter the link threshold during 3D linkage:');
    session_info.link_threshold = input('');
    
    % make a new directory for the stack in proofreading_sessions/ if necessary
    cd(proofreading_sessions_dir);
    if(exist(['./', session_info.stack_name, '/'], 'dir')~=7)
      mkdir2(['./', session_info.stack_name, '/']);
    end;
    cd(['./', session_info.stack_name, '/']);

    % make a new directory for the region in proofreading_sessions/ if necessary
    if(exist(['./', session_info.region_name, '/'], 'dir')~=7)
      mkdir2(['./', session_info.region_name, '/']);
    end;
    cd(['./', session_info.region_name, '/']);

    % make a new directory for the reconstruction
    % proofreading_session/stack/regon/ if necessary 
    if(exist(['./', session_info.reconstruction_name, '/'], 'dir')~=7)
      mkdir2(['./', session_info.reconstruction_name, '/']);
    end;
    cd(['./', session_info.reconstruction_name, '/']);
    
    % make a new directory for the new session
    mkdir2(['./', session_info.name, '/']);
    system(['chmod a+rw ', ['./', session_info.name, '/']]);
    cd(proofreading_sessions_dir);
    
    session_path_name = [proofreading_sessions_dir, session_info.stack_name, '/', ...
      session_info.region_name, '/', session_info.reconstruction_name, '/', ...
      session_info.name, '/'];
    
    % save for later access
    save2([session_path_name, session_info.name, '.session_info.mat'], '-STRUCT', 'session_info');
    
  case 4 % Continue previous session
    session_info = prev_session_info;
    session_path_name = [proofreading_sessions_dir, session_info.stack_name, '/', ...
      session_info.region_name, '/', session_info.reconstruction_name, '/', ...
      session_info.name, '/'];
    
end

fprintf('Loading the reconstruction data ..');
cd(proofreading_sessions_dir);
reconstruction_info = load2([tobe_proofread_dir, session_info.stack_name, '/', ...
  session_info.region_name, '/', session_info.reconstruction_name, '/', ...
  session_info.reconstruction_name, '.to_be_proofread.mat']);
load2(reconstruction_info.al_file);
load2(reconstruction_info.cat_file);
load2(reconstruction_info.superpixel_2_seg_map_file);
load2(reconstruction_info.links_3D_file);
reconstruction_config = reconstruction_info.config;
clear config
if(isempty(al) || isempty(cat) || isempty(superpixel_2_seg_map) || isempty(links_3D))
  error('em_reconstruction:start_proofread_session:no_recon_file', 'Could not load reconstruction, exiting.\n');
end;
fprintf('done\n\n');

al_temp = al;
al = [];
for i = 1:length(cat)
  al{i} = al_temp{i};
end;
clear al_temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step II: Launch Proofreading GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Step IV: Launching GUI\n');

global debug proof link_threshold

debug = []; proof = [];
link_threshold = session_info.link_threshold;

cd(session_path_name);
fprintf('--------- Proofreader -------------\n');
hGui = gui();
waitfor(hGui);
fprintf('--------- Exiting proofreader -----\n\n');

fprintf('Session completed\n\n');
cd(proofreading_sessions_dir);
save2('_start_proofread_session.log.mat', '-STRUCT', 'session_info');
return
