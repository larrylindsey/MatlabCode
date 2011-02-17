function [x y] = reconstructDomainBounds(instr)

if isfield(instr, 'tr')
    trinfo = instr;
elseif isfield(instr, 'section')
    trinfo = collectReconstructTransforms(instr);
else
    error('Unrecognized input struct');
end

x = zeros(2, numel(trinfo));
y = zeros(2, numel(trinfo));

for i_t = 1:numel(trinfo)
    ctrinfo = trinfo(i_t);
    
    u = ctrinfo.u(:);
    v = ctrinfo.v(:);
    
    U = cat(2, cat(1, u, u(end:-1:1)), cat(1, v, v));
    X = ctrinfo.tr.forward_fcn(U, ctrinfo.tr.tdata);
    
    x(1, i_t) = min(X(:,1));
    x(2, i_t) = max(X(:,1));
    y(1, i_t) = min(X(:,2));
    y(2, i_t) = max(X(:,2));
end
