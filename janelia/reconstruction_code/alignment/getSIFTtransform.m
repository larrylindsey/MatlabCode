function [D, err, final_match_points] = ...
  getSIFTtransform(frames,matches,inds,adj,I,options, match_points_planB)
% [D,err]=getSIFTtransform(frames,matches,inds,alpha,I,options);
% Calculates affine transform based on matches using
% SIFT landmarks and RANSAC.
% INPUTS:
%   frames  - 1xK cell-array of 4xN SIFT landmarks
%   matches - KxK cell-array of 2xM match indices into frames
%   inds    - back-reference array of indices
%   alpha   - KxK suggested adjacency matrix
%   I       0 1xK cell array of images
%   options - structure of parameters
%       .w - fraction of best matches under current model
%            to be considered inliers [0.25]
%       .p - matches with deviation under current model within
%            this factor of current inliers' mean deviation
%            to be added to inliers [2]
%       .nthr - allowed discrepancy between consistent matches [10]
%       .Nthr - minimal # of matches for evaluation of consistency [10]
%       .athr - fraction of consistent matches to declare adjacent [0.5]
%       .cthr - minimal number of RANSAC iterations [3]
%       .rank_perturb_noise - amount of noise added to vote count of
%                             consistency. This induces a certain amount of
%                             randomness to the ranking, and thus the
%                             sampling [7].
%   match_points_planB  point correspondences to be used in case SIFT
%                         feature point matching doesn't provide adequate
%                         number of matches. 
% OUTPUTS:
%   D       - 2xK cell array of found affine transforms
%             D={a{1..K},b{1..K}} such that x{k}->a{k}*x{k}+b{k}
%   err     - mean error in last RANSAC iteration
%   final_match_points    - inlier correspondence points
%
% USAGE:
% (0) make I into cell-array {plane,tiles} of images.
%
% (1) [frames,descriptors,matches]=getSIFTmatches(images,adj);
% (2) [frames1,matches1,ind1]=getSIFTrow(frames,matches);
% (3) D=getSIFTtransform(frames1,matches1,ind1);
%
% Composed based on A. Vedaldi package by Y.Mishchenko Chklovskii Lab 2008
if(nargin<3 || isempty(inds))
  inds=1:length(frames);
end
if(nargin<4)
  adj=[];
end
if(nargin<5)
  I={};
end
if(nargin<6)
  options=[];
end
if(nargin<7)
  match_points_planB=[];
end
if(isfield(options,'w'))
  w=options.w;
else
  w=0.25;
end
if(isfield(options,'p'))
  p=options.p;
else
  p=2;
end
if(isfield(options,'nthr'))
  nthr=options.nthr;
else
  nthr=10;
end
if(isfield(options,'Nthr'))
  Nthr=options.Nthr;
else
  Nthr=10;
end
if(isfield(options,'athr'))
  athr=options.athr;
else
  athr=0.3;
end
if(isfield(options,'cthr'))
  cthr=options.cthr;
else
  cthr=5;
end
if(isfield(options,'rank_perturb_noise'))
  rank_perturb_noise = options.rank_perturb_noise;
else
  rank_perturb_noise = 5;
end

if(isfield(options,'SIFT') && isfield(options.SIFT,'RANSAC') ...
    && isfield(options.SIFT.RANSAC,'consistency') ...
      && ~isempty(options.SIFT.RANSAC.consistency))
  athr=options.SIFT.RANSAC.consistency;
end
%make this for verbouse
v=0;
if(v==1)
  close all;
end

ndim=2;                     %dimensionality: !!! ndim==2 !!!
ntiles=length(frames);      %number of images to stitch

%CONSTRUCT PLANES REFERENCE
planes=zeros(1,length(frames));
for plane1=1:size(inds,1)
  for tile1=1:size(inds,2)
    if(~iscell(inds))
      % default: one patch per tile
      cnt1=inds(plane1,tile1);
    else
      % multiple patches per tile, e.g., in presence of folds
      cnt1=inds{plane1,tile1};
    end
    if(cnt1~=0)
      planes(cnt1)=plane1;
    end
  end
end

%CONSTRUCT ARRAYS OF POINTS for each pair-of-tiles m x n
x0=cell(ntiles,ntiles); y0=cell(ntiles,ntiles);
for n=1:ntiles
  for m=1:ntiles
    %need to check empty or get error
    if(isempty(matches{m,n}))
      x0{m,n}=zeros(2,0);
      y0{m,n}=zeros(2,0);
    else
      x0{m,n}=frames{m}(1:ndim,matches{m,n}(1,:));
      y0{m,n}=frames{n}(1:ndim,matches{m,n}(2,:));
    end
  end
end


% RANSAC -- FIGHT OUTLIERS -- NULL STARTUP
%null-initialization of transforms
a=cell(1,ntiles);
for n=1:ntiles
  a{n}=diag(ones(1,ndim));
end
b=cell(1,ntiles);
for n=1:ntiles
  b{n}=zeros(ndim,1);
end
inliers=cell(ntiles,ntiles);
alpha=zeros(ntiles,ntiles);
for n=1:ntiles
  for m=n+1:ntiles
    if(~isempty(planes) && abs(planes(n)-planes(m))>1)
      alpha(m,n)=0;
      alpha(n,m)=0;
      continue;
    end

    if(length(x0{m,n})<Nthr)
      alpha(m,n)=0;
      alpha(n,m)=0;
      continue;
    end

%     d=dist(y0{m,n}-x0{m,n});
    d = distance_rotation_invariant(y0{m,n}, x0{m,n});

    %fraction of consistent matches
    nfix=sum(d(:)<nthr,1);
    alpha(m,n)=sqrt(nfix/numel(d));
    alpha(n,m)=sqrt(nfix/numel(d));
    
    % bias random sampling towards consistent matches
    votes = sum(d<nthr,2);
    votes = rank_perturb_noise*(rand(size(x0{m,n},2),1)-0.5) + votes;
    [votes_sorted, sort_id] = sort(votes, 1, 'descend');
    x0{m,n} = x0{m,n}(:, sort_id);
    y0{m,n} = y0{m,n}(:, sort_id);
    x0{n,m} = x0{n,m}(:, sort_id);
    y0{n,m} = y0{n,m}(:, sort_id);
  end
end
clear d d1 idx idxa idxb
if(v==1)
  fprintf('Found consistencies:\n');
  display(alpha);
end
%estimate adjacency - make sure there is sufficient number of points
alpha=double(alpha>athr);

alpha(adj==0)=0;
%drop disjoint subset
id=find(sum(alpha,2)==0);
if(~isempty(id) && isempty(adj))
  fprintf('WARNING: found disjoint segment for SIFT in ');
  fprintf('adjacencies:\n');
  display(alpha);
  fprintf('Ids: '); fprintf('%i,',id); fprintf('\n');
  fprintf('Quitting...\n');
  D={}; err=Inf; final_match_points = [];
  return;
end

%random inliers initialization
for n=1:ntiles
  for m=n+1:ntiles
    if(isfield(options, 'n_initial_samples'))
      if(size(x0{n,m},2)<options.n_initial_samples)
        inliers{n,m} = [];
      else
        flags = zeros(1,size(x0{n,m},2));
        flags(1:options.n_initial_samples) = 1;
        inliers{n,m} = flags>0;
      end
    else
      % w fraction of matches considered as inliers
      inliers{n,m}=false(1, size(x0{n,m},2));
      inliers{n,m}(1:floor(size(x0{n,m},2)*w)) = true;
      if(size(x0{n,m},2)>6)
        % ensure that there are at least 6 inliers
        inliers{n,m}(1:6) = true;
      end
    end
    inliers{m,n}=inliers{n,m};
  end
end

%%%
% The number of SIFT detections may be inadequate when the area of overlap
% is small. For such cases use a correlation based correspondence
% detection.
% Check to make sure that each tile is atleast connected to minimum number
% of tiles. If number of tiles is greater than 2 then try to connect each
% tile to atleast two other tiles.
%%%
if(~isempty(match_points_planB))
  for tile_id_1 = 1:size(alpha,1)
    if((size(alpha,2)>2 && nnz(alpha(tile_id_1,:))<=1) || ...
        (size(alpha,2)<=2 && nnz(alpha(tile_id_1,:))<=0))
      for tile_id_2 = 1:size(alpha,2)
        if(alpha(tile_id_1, tile_id_2)>0)
          continue;
        end
        if(~isempty(match_points_planB(tile_id_1, tile_id_2).points_1))
          x0{tile_id_1, tile_id_2} = match_points_planB(tile_id_1, tile_id_2).points_1;
          y0{tile_id_1, tile_id_2} = match_points_planB(tile_id_1, tile_id_2).points_2;
          inliers{tile_id_1, tile_id_2} = true(1, size(x0{tile_id_1, tile_id_2},2));
          alpha(tile_id_1, tile_id_2) = 2;
        end
      end
    end
  end
end

if(v==1)
  fprintf('Using adjacencies:\n');
  display(alpha);
end


%some auxiliary stuff
err=Inf;
derr=Inf;
h=cell(ntiles,ntiles);

cnt=0;
while((derr>0 && err>0.1 && err/derr<100) || isinf(err) || cnt<=cthr)
  cnt=cnt+1;
  %GET TOGETHER THE MATCHING PAIRS
  x1=cell(ntiles,ntiles);
  y1=cell(ntiles,ntiles);
  for n=1:ntiles
    for m=1:ntiles
      x1{n,m}=x0{n,m}(:,inliers{n,m});
      y1{n,m}=y0{n,m}(:,inliers{n,m});
    end
  end

  %FIXED TILE - tile with the largest number of neighbors
  a=sum(alpha,2);
  [a,nfix]=max(a);

  %MODEL ESTIMATION FROM MATCHES IS HERE
%   nfix=2;
  [a,b]=getSIFTmodel(x1,y1,alpha,[],nfix);
%   [a,b]=getSIFTmodel_reg(x1,y1,alpha,[],nfix, 1);

%   for i = 1:length(a)
%     fprintf('### %g ###\n', 4*a{i}(1,2)*a{i}(2,1) + (a{i}(1,1)-a{i}(2,2)).^2);
%   end

  %++++ RANSAC RELATED  ++++
  %GET INLIERS FOR RANSAC ITERATION
  ds1=zeros(ntiles,ntiles);   %distribution of errors over tile-pairs
  for n=1:ntiles
    for m=n+1:ntiles
      if(alpha(n,m)~=1)
        continue;
      end
      %model transformed points
      x1=a{n}*x0{n,m}+repmat(b{n},[1 size(x0{n,m},2)]);
      y1=a{m}*y0{n,m}+repmat(b{m},[1 size(y0{n,m},2)]);

      %model matches discrepancies
      ds=sum((x1-y1).^2,1);
      ds1(n,m)=mean(ds(inliers{n,m})); ds1(m,n)=ds1(n,m);

      %choose inliers
      ds2=sort(ds);
      %(1) lower quartile
      ds_thr=ds2(ceil(length(ds2)*w));
      if(cnt>cthr)  %(2) and STD
        ds_thr=max(ds_thr,p*std([-ds1(n,m),ds1(n,m)]));
      end

      %these are our inliers
      inliers{n,m}=ds<=ds_thr;
      inliers{m,n}=inliers{n,m};

      %verbose
      if(v==1 && alpha(n,m)>0 && ~isempty(I))
        if(isempty(h{n,m}))
          h{n,m}=figure;
        else
          figure(h{n,m}); clf(h{n,m});
        end
        plotmatches(I{n},I{m},frames{n}(1:2,:),frames{m}(1:2,:),...
          matches{n,m}(:,inliers{n,m}));
        title(sprintf('%i<->%i(adj%i) N=%i',n,m,alpha(n,m),...
          sum(inliers{n,m})));
        pause(0.5)
      end
    end
  end

  %some cosmetics for verbose
  if(cnt==1 && v==1)
    msgbox('Arrange windows as you wish and press any key.');
    pause;
  end

  %mean inlier error, for stopping
  dst=ds1(alpha>0);
  derr=err-sqrt(mean(dst));
  err=sqrt(mean(dst));

  if(v==1)
    fprintf('RANSAC Error %4.4g/%4.4g\n',err,derr);
  end
end

%OUTPUT
D=[a,b];

final_match_points = [];
for n = 1:ntiles
  for m = 1:ntiles
    final_match_points(n,m).points_1 = x0{n,m}(:, inliers{n,m});
    final_match_points(n,m).points_2 = y0{n,m}(:, inliers{n,m});
  end
end
