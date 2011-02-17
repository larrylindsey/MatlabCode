if(isunix)
  curr_dir = pwd2;
  cd('../../lib/fftw-3.2/');
  system('./configure CFLAGS="-fPIC" --prefix=$PWD --exec-prefix=$PWD');
  system('make');
  system('make install');
  cd(curr_dir);
end
