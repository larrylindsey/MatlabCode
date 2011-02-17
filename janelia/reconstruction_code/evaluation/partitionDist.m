function [a_not_b,b_not_a] = partitionDist( a, b )
% Rand coefficient
% treats all labels equally, including zero
n = length(a(:));
if n~=length(b(:)), error('Partitions should be of the same size'); end

a_not_b = countSplitBy( a, b ) / (n*(n+1)/2);
b_not_a = countSplitBy( b, a ) / (n*(n+1)/2);

function N = countSplitBy( a, b )
    N=0;
    [asort,ix] =sort(a(:));
    marks=[0; find(diff(asort)); length(asort)];
    bSa = b(ix);
    for ia=1:(length(marks)-1)
        arr=bSa( (marks(ia)+1) : marks(ia+1) );
        cl = classCounts( arr );
        N = N +( length(arr)^2 - sum(cl.^2) )/2;
    end

function cl=classCounts( arr )
    arr = sort( arr );
    fda=[find(diff(arr)~=0)' length(arr)];
    cl=[fda(1) diff(fda)]; %class counts within b_ia

%{
a_not_b=0;
[asort,ix] =sort(a(:));
marks=[0; find(diff(asort)); length(asort)];
bSa = b(ix);
for ia=1:(length(marks)-1)
    arr=bSa( (marks(ia)+1) : marks(ia+1) );
    cl = classCounts( arr );
    a_not_b = a_not_b +( length(arr)^2 - sum(cl.^2) )/2;
end
a_not_b = a_not_b / n / (n+1) / 2;
%}