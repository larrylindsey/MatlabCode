function to_be_proofread_dir = get_to_be_proofread_dir(config)
% image_dir = get_to_be_proofread_dir(config)
%
% Get the stack image directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

stack_config = config.stack;
reconstruction_config = config.reconstruction;
to_be_proofread_config = config.to_be_proofread;

prev_dir = pwd2;
cd(to_be_proofread_config.root_dir);
if(exist(stack_config.name, 'dir')~=7)
  mkdir2(stack_config.name)
  if(isunix)
    system(['chmod -R a+rw ', stack_config.name]); 
  end
end;
cd(stack_config.name);
if(exist(config.region.dir, 'dir')~=7)
  mkdir2(config.region.dir)
  if(isunix)
    system(['chmod -R a+rw ', config.region.dir]); 
  end
end;
cd(config.region.dir);
if(exist(reconstruction_config.name, 'dir')~=7)
  mkdir2(reconstruction_config.name)
  if(isunix)
    system(['chmod -R a+rw ', reconstruction_config.name]); 
  end
end;
cd(prev_dir);

to_be_proofread_dir = [to_be_proofread_config.root_dir, stack_config.name, '/', ...
  config.region.dir, reconstruction_config.name, '/'];

return
end