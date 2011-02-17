%
% prepare all pb's for combinations with ... whatever
%
%% models
load 'modelOEOBri-HRP1.mat'; % model 5*2, evalOEOBri 5*2, code='~/segbench/hrp1/runTrainOEOBri.m'
% models are built for uint8
[nsigma,nelong]=size(model);

%% cube - please replace with your own; image should be im(1:1000,1:1000,1:10)
% load('training_data/HRP_Alex_volume.mat', 'im');
% for i = 1:size(im, 3)
%   im(:,:,i) = histeq(im(:,:,i));
% end;

load('../HRP_evaluation/datasets/images.HRP1cube1_Alex.mat', 'images');
for i = 1:length(images)
  im(:,:,i) = im2double(images{i});
end;

%% load and apply model
rFitParab = 2.5; %const
tic;
for isigma=1:nsigma
    for is=1:size(im,3) %slices
        pb{isigma}(:,:,is) = max(applyPBmodelOEOBri( model(isigma,1), im(:,:,is), rFitParab ), [],3);
    end
end
toc % 30min
%%
%isigma=5;
%view3(im,pb{isigma}.*(0==watershed(pb{isigma},4)))
