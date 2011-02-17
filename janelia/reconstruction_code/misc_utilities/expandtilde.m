function fullpath = expandtilde(somepath)
% FULLPATH = EXPANDTILDE(SOMEPATH): expand the home directory marker "~"
% in a file path to the full path.
% Although Matlab functions (such as fopen) do this automatically, .mex
% files written in C/C++ often don't offer this convenience. When using a
% function written in C, use the expandtilde() function on your filenames.

    if strcmp(somepath(1:2), '~/'),
        fullpath = [getenv('HOME') somepath(2:end)];
    else
        fullpath = somepath;
    end
end