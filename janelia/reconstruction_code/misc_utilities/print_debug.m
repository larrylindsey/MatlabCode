function print_debug(format_string, varargin)
% print_debug(format_string, varargin)
% Prints in debug mode to standard output, otherwise silent

global config_global
if(config_global.DEBUG)
  fprintf(format_string, varargin{:});
end

return
end
