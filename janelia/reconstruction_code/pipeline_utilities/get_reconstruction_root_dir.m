function reconstruction_root_dir = get_reconstruction_root_dir()
% get_reconstruction_root_dir()
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  02092009  init code
%

prev_dir = pwd2();
curr_dir = prev_dir;
found_reconstruction_dir = false;
while((strcmp(curr_dir, '/')~=1) && ~found_reconstruction_dir)
  cd('..');
  curr_dir = pwd2();
  if(exist('reconstructions', 'dir')==7)
    found_reconstruction_dir = true;
  end
end
if(~found_reconstruction_dir)
  error('Not in correct directory. Please go within reconstructions/scripts/ and run again.');
end
reconstruction_root_dir = [pwd2, '/'];
addpath([reconstruction_root_dir, 'scripts/']);
cd(prev_dir);
return
end
