clear all
close all

temp_file_name = '/groups/chklovskii/home/vitaladevunis/temp/a.raw';

a = int8(rand(100, 1)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read int8 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read int8 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = int8(rand(100, 210)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read int8 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read int8 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = int8(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read int8 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read int8 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = uint8(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read uint8 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read uint8 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = int16(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read int16 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read int16 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = uint16(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read uint16 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read uint16 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = int32(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read int32 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read int32 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = uint32(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read uint32 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read uint32 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = int64(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read int64 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read int64 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = uint64(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read uint64 array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read uint64 array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = single(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read single array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read single array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);

a = double(rand(100, 210, 17)*997);
fwrite_raw_array_mex(temp_file_name, a);
[b, err] = fread_raw_array_mex(temp_file_name);
if(err==0 && nnz(a~=b)==0)
  fprintf(['Successfully wrote and read double array of size [', ...
    num2str(size(a)), ']\n']);
else
  error(['Unsuccessful test for and read double array of size [', ...
    num2str(size(a)), ']. Error code: ', num2str(err), '\n']);
end
delete(temp_file_name);



