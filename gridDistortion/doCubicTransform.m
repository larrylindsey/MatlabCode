function [U A trans] = doCubicTransform(X, trans, direction)
% [U A trans] = doCubicTransform(X, trans, direction)
% Performs an order three bivariate transform on the input data X
%
% X - [n 2] - the first column represents data in the first coordinate (x),
%       the second column, data in the second coordinate (y).
% trans - struct - the transform struct to be applied to X
%         FIELDS
%         trans.type - if textual, indicates either 'taylor' or 'legendre'
%                      type bases.  Otherwise, should hold a handle to a
%                      function that will return a matrix in a similar
%                      fashion to taylorMatrix or legendreMatrix
%         trans.T - [m 2] - the transform matrix, as according to the
%                           relevant basis.
%         trans.Tinv - [m 2] - if non-empty, should contain the inverse
%                              transform matrix.  This is most often
%                              approximate.
%         trans.order - [1] - this indicates the maximal order of the
%                             transform. This method is only compatible up
%                             to order 3.
%         trans.data - struct - Generally, this field is unspecified.  This
%                               method will attempt to find fields n, u,
%                               and v in the event that a reverse transform
%                               is attempted, and Tinv is empty.  This
%                               information will be used to compute an
%                               approximate inverse transform matrix to T.
% direction - logical - false for forward, true for reverse.
%
% U - [n 2] - the result of applying the transform to X
% A - [m 2] - the basis matrix, as returned by legendreMatrix,
%             taylorMatrix, etc.
% trans - struct - this function passes back the transform struct argument
%                  to allow the possibility of preserving the inverse
%                  matrix, if desired.

% Note:  The type field may in fact contain any textual information for
% which there is a correspoding Matrix function, ie
% 'legendre' -> @legendreMatrix
% 'taylor' -> @taylorMatrix
% 'foobar' -> @foobarMatrix

if isfield(trans, 'tdata');
    trans = trans.tdata;
end

if isfield(trans, 'order')
    order = trans.order;
else
    order = 3;
end

tsize = (order + 1) * (order + 2) / 2;

if nargin < 3
    direction = false;
end

if ischar(trans.type)
    type = str2func([trans.type 'Mat']);
else
    type = trans.type;
end

A = type(X(:,1), X(:,2), order);

if direction
    T = trans.Tinv;
    if isempty(T)
        T = invertTransform(trans.T, type, trans.data.n, ...
            trans.data.u);
        trans.Tinv = T;
    end
else
    T = trans.T;
end

T = T(1:tsize,:);

U = A * T;

