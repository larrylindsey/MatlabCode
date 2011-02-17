function em_reconstruct(config_param, module_ids, varargin)
% em_reconstruct(config_file_name, module_ids, varargin) executes modules
% specified in module_id list for the config parameters specified by 
% config_param.
% Inputs
%   config_param        Two options
%                       (1) A '.m' or '.mat' file name. If '.m' then must be a
%                         function returning config. If '.mat', then
%                         variable config is loaded from it.
%                       (2) A struct specifying the reconstruction
%                       parameters through two sub-structs:
%                       config_param.config and config_param.config_global.
%   module_ids          List of modules to be executed. See pipeline,
%                         pipeline_serial_section, pipeline_block_face for
%                         a list of valid module ids.
%   optional parameter/value list:
%     case_ids            List of sections to be processed [empty->do all] 
%     is_verbose          Whether to print messages [true]
%     is_verbose_figures  Whether to display figures for intermediate
%                           results [false]
%     is_stand_alone      Whether to produce stand-alone batch scripts
%                           [false]
%     job_name            Name to be given to the job.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  04212009    init. code
%

global config_global

if(isdeployed)
  module_ids = eval(module_ids);
end

% Parse the input arguments
case_ids = [];
is_verbose = true;
is_verbose_figures = false;
is_stand_alone = false;
job_name = '';
notify_email = [];
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'case_ids'
      case_ids = varargin{i+1};
      if(isdeployed)
        case_ids = eval(case_ids);
      end
    case 'is_verbose'
      is_verbose = varargin{i+1};
    case 'is_verbose_figures'
      is_verbose_figures = varargin{i+1};
    case 'is_stand_alone'
      is_stand_alone = varargin{i+1};
    case 'job_name'
      job_name = varargin{i+1};
    case 'notify_email'
      notify_email = varargin{i+1};
    otherwise
      error('Option not understood.');
  end
end

if(is_verbose)
  fprintf('START: em_reconstruct\n');
end

% Construct the config structure
if(isa(config_param, 'char'))
  if(strcmp(config_param(end-3:end), '.mat')~=1)
    config = get_basic_config();
    
    config_global.job.is_stand_alone = is_stand_alone;
    if(~isempty(job_name))
      config.job.name = job_name;
    end
    
    h_config_f = str2func(config_param);
    config = h_config_f(config, case_ids, is_verbose, is_verbose_figures);
  else
    load_config = load(config_param);
    config = load_config.config;
    config_global = load_config.config_global;
    if(~isempty(case_ids))
      config.stack.case_ids = case_ids;
    end
  end
else
  config = config_param.config;
  config_global = config_param.config_global;
  if(~isempty(case_ids))
    config.stack.case_ids = case_ids;
  end
end

config.is_verbose = is_verbose;
config.is_verbose_figures = is_verbose_figures;

% Call the EM reconstruction pipeline.
pipeline(module_ids, config);

if(~isempty(notify_email))
  notify_email_subject = ['em_reconstruct(''', config_file_name, ''', ', ...
    '[', num2str(module_ids, '%d '), ']) has executed.'];
  notify_email_text = ['em_reconstruct(''', config_file_name, ''', ', ...
    '[', num2str(module_ids, '%d '), ']) has executed.'];
  send_email_python(notify_email.gmailUser, notify_email.gmailPassword, ...
    notify_email.recipient, notify_email_subject, notify_email_text);
end

if(is_verbose)
  fprintf('STOP: em_reconstruct\n');
end

return;
end
