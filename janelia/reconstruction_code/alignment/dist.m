function D=dist(X)
%D=dist(X)
%Partial implementation of dist-function computing all
%pairwise distances for vectors in ndim x L matrix of 
%vectors X

v=sum(X.^2,1);
D=repmat(v,[size(X,2) 1]);
D=sqrt(D+D'-2*X'*X);
