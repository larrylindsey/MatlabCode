function disp_debug(varargin)
% disp_debug(format_string, varargin)
% Prints in debug mode to standard output, otherwise silent

global config_global
if(config_global.DEBUG)
  disp(varargin{:});
end

return
end
