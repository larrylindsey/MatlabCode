function trans = transformStruct(type, T, Tinv, u, v)

if nargin < 5
    v = [];
    if nargin < 4
        u = [];
        if nargin < 3
            Tinv = [];
        end
    end
end

trans.type = type;
trans.T = T;

trans.Tinv = Tinv;
   
TOrderTest = sum(abs(T), 2);
i_order = find(TOrderTest > 0, 1, 'last') - 1;

order = -1;

while i_order > 0
    i_order = i_order - order - 2;
    order = order + 1;
end

trans.order = order;
trans.data.n = 32;

trans.data.u = u;
trans.data.v = v;

trans.iDim = 2;
trans.oDim = size(T, 2);
