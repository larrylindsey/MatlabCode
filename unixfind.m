function [o s] = unixfind(varargin)
% function [o s] = unixfind([arg1, [arg2 ... ] ]);
% A function to simplify calling the unix find command.
% Passes each argument directly to find, in order
%
% o - a cell array containing the output of unix find
% s - the integer return status

o = {};
cmd = 'find';

for iarg = 1:numel(varargin)
    cmd = sprintf('%s %s', cmd, varargin{iarg});
end

[s r] = unix(cmd);

[o{1} tokr] = strtok(r, sprintf('\n'));
while ~isempty(strtrim(deblank(tokr)))
    [o{end + 1} tokr] = strtok(tokr, sprintf('\n')); %#ok<AGROW,STTOK>
end

o = strtrim(deblank(o));
