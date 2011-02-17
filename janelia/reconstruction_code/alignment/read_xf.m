function xf = read_xf(xf_file_name)
% xf = xf_2_tform(xf_file_name)
% Read in .xf file.
%
% Inputs:
%   xf_file_name  name with full path of .xf file.
% Output:
%   xf            [N x 6] matrix where N = #slices
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  05082008  init code
%

fin = fopen(xf_file_name, 'rt');
xf = [];
while(feof(fin)~=1)
  x = fscanf(fin, '%f', [1 6]);
  if(~isempty(x))
    xf(end+1,:) = x;
  end;
end

return
end
