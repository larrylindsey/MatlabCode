function [U A] = doCubicTransform2(X, T, varargin)

select = false;
type = 0; %taylor, as opposed to legendre

while numel(varargin) > 0
    arg = varargin{1};
    varargin = {varargin{2:end}};
    if isa(arg, 'logical')
        select = arg;
    elseif ischar(arg)
        if strcmp(arg, 'taylor')
            type = 0;
        elseif strcmp(arg, 'legendre')
            type = 1;
        else
            error('Unrecognized argument %s', arg);
        end
    else
        disp(arg);
        error('Unrecognized argument');
    end        
end


%X in P by 2
%T in 6 by 2 or 12 by 2

if ~isnumeric(X)
    if nargin > 1
        if isa(X, 'logical') && ischar(T)
            a = X;
            b = strcmp(T, 'legendre');
        elseif isa(T, 'logical') && ischar(X)
            a = T;
            b = strcmp(X, 'legendre');
        else
            disp(X);
            disp(T);
            error('Unusable argument list');
        end
    else
        if isa(X, 'logical')
            a = X;
            b = false;
        else
            disp(X);
            error('Unusable argument');
        end
    end


    if a
        select = 'true';
    else
        select = 'false';
    end           
    
    if b
        basis = 'legendre';
    else
        basis = 'taylor';
    end
    
    U = eval(['@(X, T) doCubicTransform2(X, T, ' select ', '''...
         basis ''');']);
else    
    %addIndex = select * 6;

    if isstruct(T)
        T = T.tdata;
    end

    if select
        T = T(11:20, :);
    else
        T = T(1:10, :);
    end

    

    switch type
        case 0
            A = taylorMatrix(X(:,1), X(:,2));
        case 1
            A = legendreMatrix(X(:,1), X(:,2));
    end
    
    U = A * T;
    
end

end
