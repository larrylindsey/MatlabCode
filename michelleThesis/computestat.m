function [stat, e] = computestat(sname, cname, str, type, minArea, cal)
% sname - stat name
% cname - class name
% str - the struct
% type - 'a' or 'm' for annotation or mito

% population
pop = str.(cname).(type);
pop(pop < minArea) = [];

if nargin > 5
    pop = pop * cal;
end

if nargout > 1
    e = stderr(pop);
end

switch sname
    case 'n'
        stat = numel(pop);
    case 'fraction'
        % non capillary pixels
        noncappix = str.totpix - sum(str.capillary.a);
        stat = sum(pop) / noncappix;
    case 'avgsize'
        stat = mean(pop);
    case 'nnorm'
        % non capillary pixels
        noncappix = str.totpix - sum(str.capillary.a);
        stat = numel(pop) / noncappix;
    otherwise
        error('don''t know how to compute stat %s', sname);
end

end