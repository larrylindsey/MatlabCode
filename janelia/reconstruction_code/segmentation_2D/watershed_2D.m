function watershed_2D(boundary_map_file_name, watershed_file_name, conn)

if(isdeployed)
  if(nargin<=2)
    conn = 4;
  else
    conn = eval(conn);
  end
else
  if(nargin<=2)
    conn = 4;
  end
end

b = imread(boundary_map_file_name);
ws = watershed(b, conn);
if(max(ws(:))>65536)
  error('Maximum watershed segment label exceeds 65536. Will result in error when saving');
end

imwrite(uint16(ws), watershed_file_name);

return
end
