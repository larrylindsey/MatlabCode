function trans = convertTransMat2Struct(T, type, u, v)

trans.type = type;
trans.T = T([10 8:9 5:7 1:4],:);

if size(T, 1) > 10
    trans.Tinv = T(10 + [10 8:9 5:7 1:4],:);
else
    trans.Tinv = [];
end
   
trans.order = 3;
trans.data.n = 32;

if nargin > 2
    trans.data.u = u;
else
    trans.data.u = [];
end

if nargin > 3
    trans.data.v = v;
else
    trans.data.v = [];
end


trans.iDim = 2;
trans.oDim = 2;
