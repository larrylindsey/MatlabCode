function trans = populateTransInverse(trans)

if isempty(trans.T)
    if ~isempty(trans.Tinv)
        trans = invertTransStruct(...
            populateTransInverse(invertTransStruct(trans)));
    end
end

T = trans.T;
% if any(T(4:end,:) ~= 0)
    if isfield(trans, 'data')
        if isfield(trans.data, 'u')
            u = trans.data.u;
        else
            u = [];            
        end
        if isfield(trans.data, 'v')
            v = trans.data.v;
        else
            v = [];
        end
        if isfield(trans.data, 'n')
            n = trans.data.n;
        else
            n = []; 
        end
    else
        u = [];
        v = [];
        n = [];
    end
    support = [min([u(:) ; v(:)]) max([u(:) ; v(:)])];
    Tinv = invertTransform(T, trans.type, n, support);
% else
%     A = T(2:3,:)';
%     Ainv = inv(A);
%     c = T(1,:)';
%     
%     cinv = - Ainv * c;
%     
%     Tinv = cat(1, cinv', Ainv', zeros(size(T, 1) - 3, 2));
%     
% end

trans.Tinv = Tinv;
