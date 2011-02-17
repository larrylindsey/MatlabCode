function [a,b]=getSIFTmodel(x,y,alpha,w,n0)
% A=getSIFTmodel(x,y,alpha,w,n0)
% generates set of affine models fitting x{m}->y{n} such that x=a*y+b
% inputs: cell-array KxK of matched points 2 x N x<->y
%         cell-array KxK of point weights 1xN w
%         array KxK of pairs weights alpha
%         fixed tile index n0
%
% Composed based on A. Vedaldi package by Y.Mishchenko Chklovskii Lab 2008
ndim=2;             %number of dimensions
ntiles=size(x,1);   %number of tiles in set
%initialize if missing anything
if(nargin<3 || isempty(alpha))
  alpha=zeros(ntiles,ntiles);
  for n=1:ntiles
    for m=1:ntiles
      alpha(m,n)=1;
    end
  end
end
if(nargin<4 || isempty(w))
  w=cell(ntiles,ntiles);
  for n=1:ntiles
    for m=1:ntiles
      w{m,n}=ones(1,size(x{m,n},2));
    end
  end
end
if(nargin<5 || isempty(n0)) n0=1; end


% CONSTRUCT EQUATIONS FOR M*[a(:);b(:)]=N
% FROM QUADRATIC MINIMIZATION PROBLEM ON AFFINE
LHS=[]; RHS=[];
%for each n
for n=1:ntiles
  fprintf('%d ', n)
  if(mod(n,20)==0)
    fprintf('\n');
  end
  %fixed tile is FIXED
  if(n==n0) continue; end

  %for each d/db^n_i
  for i=1:ndim
    M=[];
    for m=1:ntiles
      A=zeros(ndim,ndim); B=zeros(ndim,1);
      if(m==n0)
        %if hit fixed tile, this is RHS
        if(isempty(x{m,n}))
          v=0;
        else
          v=alpha(m,n)*sum(w{m,n}.*x{m,n}(i,:));
        end
        RHS=[RHS;v];
      elseif(m==n)
        %if hit yourself, this is plus-part
        s=0;
        for l=1:ntiles
          if(~isempty(y{l,n}) && l~=n)
            s=s+alpha(l,n)*sum(w{l,n});
          end
        end
        B(i)=s;

        v=zeros(2,1);
        for l=1:ntiles
          if(~isempty(y{l,n}) && l~=n)
            v=v+alpha(l,n)*sum(repmat(w{l,n},[2 1]).*y{l,n},2);
          end
        end
        A(i,:)=v(:)';
        M=[M,[A(:)',B(:)']];
      else
        %otherwise is minus-part
        s=alpha(m,n)*sum(w{m,n});
        B(i)=-s;

        if(isempty(x{m,n}))
          v=zeros(ndim,1);
        else
          v=alpha(m,n)*sum(repmat(w{m,n},[2 1]).*x{m,n},2);
        end
        A(i,:)=-v(:)';
        M=[M,[A(:)',B(:)']];
      end
    end
    LHS=[LHS;M];
  end

  %for each d/da^n_ij
  for i=1:ndim
    for j=1:ndim
      M=[];
      for m=1:ntiles
        A=zeros(ndim,ndim); B=zeros(ndim,1);
        if(m==n0)
          %if hit fixed tile, this is RHS
          if(isempty(x{m,n}))
            v=0;
          else
            v=alpha(m,n)*sum(w{m,n}.*x{m,n}(i,:).*y{m,n}(j,:));
          end
          RHS=[RHS;v];
        elseif(m==n)
          %if hit yourself, this is plus LHS
          s=0;
          for l=1:ntiles
            if(~isempty(y{l,n}) && l~=n)
              s=s+alpha(l,n)*sum(w{l,n}.*y{l,n}(j,:));
            end
          end
          B(i)=s;

          v=zeros(2,1);
          for l=1:ntiles
            %check empty or error from .*
            if(~isempty(y{l,n}) && l~=n)
              v=v+alpha(l,n)*sum(repmat(w{l,n}.*y{l,n}(j,:),[2 1]).*y{l,n},2);
            end
          end
          A(i,:)=v(:)';
          M=[M,A(:)',B(:)'];
        else
          %otherwise is minus LHS
          %check empty or error from .*
          if(isempty(y{m,n}))
            s=0;
          else
            s=alpha(m,n)*sum(w{m,n}.*y{m,n}(j,:));
          end
          B(i)=-s;

          %check empty or error from .*
          if(isempty(y{m,n}))
            v=zeros(2,1);
          else
            v=alpha(m,n)*sum(repmat(w{m,n}.*y{m,n}(j,:),[2 1]).*x{m,n},2);
          end
          A(i,:)=-v(:)';
          M=[M,A(:)',B(:)'];
        end
      end
      LHS=[LHS;M];
    end
  end

end     %end for each n

%SOLVE LIN-SYSTEM
A=linsolve(LHS,RHS);

%DECONVOLVE SOLUTION INTO MATRIX FORM
a=cell(1,ntiles); b=cell(1,ntiles); K=ndim^2+ndim;
cnt=0;
for n=1:ntiles
  if(n==n0)
    %fixed tile
    a{n}=diag(ones(1,ndim)); b{n}=zeros(ndim,1);
  else
    %other tiles
    AA=A(K*cnt+1:K*(cnt+1));
    a{n}=zeros(ndim,ndim); a{n}(:)=AA(1:ndim^2);
    b{n}=zeros(ndim,1); b{n}(:)=AA(ndim^2+1:end);

    cnt=cnt+1;
  end
end
