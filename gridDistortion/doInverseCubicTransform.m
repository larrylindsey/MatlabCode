function [U A] = doInverseCubicTransform(X, trans)
% [U A] = doInverseCubicTransform(X, trans)
% This function is written for convenience when used with other functions
% of the type imtransform.  It is equivalent to a call to
%
% doCubicTransform(X, trans, true)

[U A] = doCubicTransform(X, trans, true);