function [err, err_msg] = send_done_email_python(id)
% [err, err_msg] = send_email_python(gmailUser, gmailPassword, recipient, ...
%   subject, text, attachments)
% Send email via python gmail
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(nargin<1)
  id = mfilename();
  id_f = dbstack();
else
  id_f = id;
end

global code_dir

global config_global

sys_com = ['python -c "', ...
  'import sys; ', ...
  'sys.path = sys.path + [''', code_dir, 'lib/python_email'']; ', ...
  'import py_email; ', ...
  'py_email.sendMail(', ...
  '\"', config_global.notify_email.gmailUser, '\"', ...
  ',\"', config_global.notify_email.gmailPassword, '\"', ...
  ',\"', config_global.notify_email.recipient, '\"', ...
  ',\"Done: ', id, ' \"', ...
  ',\"Done: ', id_f, ' \"', ...
  ');"'];

[err, err_msg] = system(sys_com);

return;
end
