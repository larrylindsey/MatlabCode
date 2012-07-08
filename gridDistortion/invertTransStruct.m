function trans = invertTransStruct(trans, prefun, postfun)
% Inverts a transform struct
%  trans = invertTransStruct(trans, [prefun, [postfun]])
%    trans - the transform to invert
%    prefun - if included and nonempty, a function handle for
%             pre-processing the transform. Must accept only a struct, and
%             output only a struct.
%    postfun - like prefun, for post-processing.

standardFields = {'T', 'doTrans', 'createTrans', 'matrixFun', 'param', 'order', ...
    'monomialOrder', 'ndim'};

if nargin > 1 && ~isempty(prefun)
        trans = prefun(trans);   
end
       
if isempty(trans.inv)
    trans = fitInverseTransform(trans);
end

for ii = 1:numel(standardFields)
    origField = trans.(standardFields{ii});
    trans.(standardFields{ii}) = trans.inv.(standardFields{ii});
    trans.inv.(standardFields{ii}) = origField;
end

if isfield(trans, 'data')
    origField = trans.data;
    trans.data = trans.inv.data;
    trans.inv.data = origField;
end    

if nargin > 2 && ~isempty(postfun)
    trans = postfun(trans);
end
