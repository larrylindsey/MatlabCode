function dmeshScript(index)

cd /media/External/Images/0'56D Series JPGS'/

addpath(genpath('~/code/matlab'))

load ~/code/matfiles/session_05_28_2010.mat

dmesh = contourGrid(secdoc(index:4:end), 'd01', '', int16(4), 512, 5);

save(sprintf('/home/larry/code/matfiles/dmesh_%d.mat', index), 'dmesh');