function sendmsg(varargin)

dest = getenv('msgdest');
if ~isempty(dest)
    try
        sendmail(dest, varargin{:});
    catch sendmsgerr
        fprintf('Unable to send message %s to %s\n', varargin{1}, dest);
        fprintf('Error message was %s\n', sendmsgerr.message);        
    end
end