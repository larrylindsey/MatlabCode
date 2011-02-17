function make_mex(current_dir)
% Build all mex files in subdirectories
%

if(nargin==0)
  current_dir = [pwd2(), '/'];
end

fprintf([current_dir, '\n']);

if(exist('./make_mex_here.m', 'file')==2)
  make_mex_here();
end

files = dir('.');
for i = 1:length(files)
  if(strcmp(files(i).name(1), '.')==1)
    continue;
  end
  
  if(files(i).isdir)
    cd(['./', files(i).name, '/']);
    make_mex([current_dir, files(i).name, '/']);
    cd(current_dir);
  end
end

return
end
